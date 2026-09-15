from pydantic import BaseModel, Field

from app.models.tumor_types import TumorType


class BatchRequest(BaseModel):
    variants: list[str] = Field(
        min_length=1,
        description="List of variants to normalize. Each variant must be in HGVS format.",
        examples=[["NM_004119.3:c.2073T>G", "NC_000023.11:g.32389644G>A"]],
    )
    tumorType: TumorType | None = Field(
        default=None,
        description="Optional tumor type used for context-dependent evidence rules.",
        examples=["Breast Cancer"],
    )
