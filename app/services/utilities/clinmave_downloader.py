"""Download and refresh locally cached ClinMAVE data files.

This utility exists because ClinMAVE does not expose a simple bulk data archive
for all per-gene variant exports. Gene-level variant CSVs are generated through an
asynchronous task API behind the Download page. This module reproduces that
browser workflow so the repository can refresh its ClinMAVE cache without
manual clicking.

What it downloads:

- A gene list resolved from the ClinMAVE gene autocomplete endpoint.
- One per-gene CSV export for each requested gene, named ``variants.<GENE>.csv``.
- A manifest recording byte counts and success or error status for each gene.

The current per-gene download flow is:

1. Query ``/clinmave/api/select/genes`` to resolve valid gene symbols.
2. Start an export job via ``/clinmave/api/download?geneName=<GENE>``.
3. Poll ``/clinmave/api/download/progress`` until the task reaches 100%.
4. Fetch the generated CSV from ``/clinmave/api/download/file``.

This module intentionally uses only the Python standard library so it remains
usable in minimal environments and so a future refresh is not blocked on extra
package installation. If ClinMAVE changes or disappears, this file preserves the
last known working integration path and makes that breakage explicit.

Usage examples:

    python -m app.services.utilities.clinmave_downloader
    python -m app.services.utilities.clinmave_downloader --gene BRAF --gene TP53
    python -m app.services.utilities.clinmave_downloader --skip-existing
"""

from __future__ import annotations

import argparse
import csv
import json
import ssl
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


CLINMAVE_BASE_URL = "https://ngdc.cncb.ac.cn/clinmave"
CLINMAVE_API_BASE_URL = f"{CLINMAVE_BASE_URL}/api"
DEFAULT_TIMEOUT_SECONDS = 120
GENE_POLL_INTERVAL_SECONDS = 1.0

SUMMARY_FILES: dict[str, str] = {}

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "clinmave"
DEFAULT_GENES_FILE = DEFAULT_OUTPUT_DIR / "genes.txt"
DEFAULT_MANIFEST_FILE = DEFAULT_OUTPUT_DIR / "download-manifest.csv"


@dataclass(slots=True)
class DownloadResult:
    gene: str
    filename: str
    bytes_written: int
    status: str


class ClinMaveClient:
    def __init__(self, timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS) -> None:
        self.timeout_seconds = timeout_seconds
        self.ssl_context = ssl.create_default_context()
        self.headers = {"User-Agent": "Mozilla/5.0"}

    def _request(self, url: str, params: dict[str, Any] | None = None) -> bytes:
        if params:
            query = urllib.parse.urlencode(params, doseq=True)
            url = f"{url}?{query}"
        request = urllib.request.Request(url, headers=self.headers)
        with urllib.request.urlopen(
            request,
            timeout=self.timeout_seconds,
            context=self.ssl_context,
        ) as response:
            return response.read()

    def get_text(self, url: str, params: dict[str, Any] | None = None) -> str:
        return self._request(url, params).decode("utf-8", "ignore")

    def get_json(self, url: str, params: dict[str, Any] | None = None) -> Any:
        return json.loads(self.get_text(url, params))

    def download_file(self, url: str, destination: Path) -> int:
        destination.write_bytes(self._request(url))
        return destination.stat().st_size

    def fetch_gene_options(self, gene_name: str = "") -> list[dict[str, Any]]:
        payload = self.get_json(
            f"{CLINMAVE_API_BASE_URL}/select/genes",
            {"geneName": gene_name},
        )
        if isinstance(payload, list):
            return payload
        if isinstance(payload, dict) and isinstance(payload.get("data"), list):
            return payload["data"]
        raise ValueError("Unexpected ClinMAVE gene list response shape")

    def start_gene_download(self, gene_name: str) -> str:
        task_id = self.get_text(
            f"{CLINMAVE_API_BASE_URL}/download",
            {"geneName": gene_name},
        ).strip().strip('"')
        if not task_id:
            raise ValueError(f"Failed to start ClinMAVE download for {gene_name}")
        return task_id

    def wait_for_download(self, task_id: str) -> None:
        for _ in range(self.timeout_seconds):
            payload = self.get_json(
                f"{CLINMAVE_API_BASE_URL}/download/progress",
                {"taskId": task_id},
            )
            progress = payload
            if isinstance(payload, dict):
                progress = payload.get("progress", payload.get("data", payload))
            if float(progress) >= 100:
                return
            time.sleep(GENE_POLL_INTERVAL_SECONDS)
        raise TimeoutError(f"Timed out waiting for ClinMAVE task {task_id}")

    def fetch_gene_csv(self, gene_name: str, destination: Path) -> int:
        task_id = self.start_gene_download(gene_name)
        self.wait_for_download(task_id)
        destination.write_bytes(
            self._request(
                f"{CLINMAVE_API_BASE_URL}/download/file",
                {"taskId": task_id},
            )
        )
        return destination.stat().st_size


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download ClinMAVE summary files and per-gene variant CSVs.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Defaults to {DEFAULT_OUTPUT_DIR}.",
    )
    parser.add_argument(
        "--gene",
        action="append",
        dest="genes",
        default=[],
        help="Gene symbol to download. Repeat to download multiple genes.",
    )
    parser.add_argument(
        "--genes-file",
        type=Path,
        default=None,
        help="Optional file with one gene symbol per line.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip files that already exist in the output directory.",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Deprecated no-op retained for CLI compatibility.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=2,
        help="Retries per gene after the initial attempt. Defaults to 2.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=DEFAULT_TIMEOUT_SECONDS,
        help=f"HTTP timeout and task wait budget per request. Defaults to {DEFAULT_TIMEOUT_SECONDS}.",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def ensure_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_gene_list(client: ClinMaveClient, args: argparse.Namespace) -> list[str]:
    genes = [gene.strip() for gene in args.genes if gene and gene.strip()]
    if args.genes_file:
        genes.extend(read_gene_file(args.genes_file))
    if not genes:
        genes = [option["value"] for option in client.fetch_gene_options("") if option.get("value")]
    unique_genes = sorted(set(genes))
    return unique_genes


def read_gene_file(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def write_gene_file(path: Path, genes: list[str]) -> None:
    path.write_text("\n".join(genes) + "\n", encoding="utf-8")


def download_summary_files(
    client: ClinMaveClient,
    output_dir: Path,
    skip_existing: bool,
) -> None:
    for filename, url in SUMMARY_FILES.items():
        destination = output_dir / filename
        if skip_existing and destination.exists():
            print(f"skip summary {filename}")
            continue
        bytes_written = client.download_file(url, destination)
        print(f"downloaded summary {filename} ({bytes_written} bytes)")


def download_gene_with_retries(
    client: ClinMaveClient,
    gene: str,
    output_dir: Path,
    skip_existing: bool,
    retries: int,
) -> DownloadResult:
    filename = f"variants.{gene}.csv"
    destination = output_dir / filename
    if skip_existing and destination.exists():
        return DownloadResult(gene=gene, filename=filename, bytes_written=destination.stat().st_size, status="ok")

    attempts = retries + 1
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            bytes_written = client.fetch_gene_csv(gene, destination)
            return DownloadResult(gene=gene, filename=filename, bytes_written=bytes_written, status="ok")
        except Exception as exc:  # pragma: no cover - network variability
            last_error = exc
            print(f"attempt {attempt}/{attempts} failed for {gene}: {exc}", file=sys.stderr)
            if destination.exists() and destination.stat().st_size == 0:
                destination.unlink()
            if attempt < attempts:
                time.sleep(GENE_POLL_INTERVAL_SECONDS)
    return DownloadResult(
        gene=gene,
        filename=filename,
        bytes_written=0,
        status=f"error:{last_error}",
    )


def write_manifest(path: Path, results: list[DownloadResult]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "filename", "bytes", "status"])
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "gene": result.gene,
                    "filename": result.filename,
                    "bytes": result.bytes_written,
                    "status": result.status,
                }
            )


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.out_dir.resolve()
    ensure_output_dir(output_dir)

    client = ClinMaveClient(timeout_seconds=args.timeout_seconds)
    download_summary_files(client, output_dir, args.skip_existing)
    if args.summary_only:
        return 0

    genes = load_gene_list(client, args)
    write_gene_file(output_dir / "genes.txt", genes)

    results: list[DownloadResult] = []
    for gene in genes:
        result = download_gene_with_retries(
            client=client,
            gene=gene,
            output_dir=output_dir,
            skip_existing=args.skip_existing,
            retries=args.retries,
        )
        results.append(result)
        print(f"{result.gene}: {result.status} ({result.bytes_written} bytes)")

    write_manifest(output_dir / "download-manifest.csv", results)
    error_count = sum(1 for result in results if result.status != "ok")
    if error_count:
        print(f"completed with {error_count} gene download errors", file=sys.stderr)
        return 1
    print(f"completed successfully: {len(results)} gene CSVs")
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entrypoint
    raise SystemExit(main())
