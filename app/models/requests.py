from pydantic import BaseModel, ConfigDict, Field

from app.models.tumor_types import TumorType


class BatchRequest(BaseModel):
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "variants": [
                    "NM_001002295.2:c.3G>A",
                    "NC_000023.11:g.32389644G>A",
                    "NC_000013.11:28027221:A:C",
                ],
                "tumorType": "Non-Small Cell Lung Cancer",
            }
        }
    )

    variants: list[str] = Field(
        min_length=1,
        description=(
            "List of variants to normalize. Each variant must be in "
            "HGVS or SPDI format on an NC_ or NM_ reference sequence."
        ),
    )
    tumorType: TumorType | None = Field(
        default=None,
        description="Optional tumor type used for context-dependent evidence rules.",
    )
