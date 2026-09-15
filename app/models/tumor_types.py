from typing import Literal, get_args


TumorType = Literal[
    "Bladder Carcinoma",
    "Breast Cancer",
    "Hodgkin Lymphoma",
    "Neuroblastoma",
    "Parathyroid Carcinoma",
    "Peripheral T-Cell Lymphoma",
    "Renal Cell Carcinoma",
    "T-Cell Acute Lymphoblastic Leukemia",
    "T-Cell Lymphoblastic Lymphoma",
    "Urothelial Carcinoma",
]

TUMOR_TYPE_VALUES = get_args(TumorType)
