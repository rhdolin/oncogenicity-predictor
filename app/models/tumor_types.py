from typing import Literal, get_args


TumorType = Literal[
    "Bladder Carcinoma",
    "Breast Cancer",
    "Granulosa Cell Tumor",
    "Hairy Cell Leukemia",
    "Hodgkin Lymphoma",
    "Neuroblastoma",
    "Non-Small Cell Lung Cancer",
    "Parathyroid Carcinoma",
    "Peripheral T-Cell Lymphoma",
    "Retinoblastoma",
    "Renal Cell Carcinoma",
    "T-Cell Acute Lymphoblastic Leukemia",
    "T-Cell Lymphoblastic Lymphoma",
    "Urothelial Carcinoma",
]

TUMOR_TYPE_VALUES = get_args(TumorType)
