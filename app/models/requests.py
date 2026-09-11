from pydantic import BaseModel, Field


class BatchRequest(BaseModel):
    variants: list[str] = Field(
        min_length=1,
        description="List of variants to normalize. Each variant must be in HGVS format.",
        examples=[["NM_004119.3:c.2073T>G", "NC_000023.11:g.32389644G>A"]],
    )
