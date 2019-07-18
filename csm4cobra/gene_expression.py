import numpy as np
import pandas as pd


TISSUES = ("ADRENAL_CORTEX", "AUTONOMIC_GANGLIA", "BILIARY_TRACT", "BONE", "BREAST",
           "BUCCAL", "CENTRAL_NERVOUS_SYSTEM", "CERVIX", "ENDOMETRIUM", "EYE", "FIBROBLAST",
           "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "KIDNEY", "LARGE_INTESTINE", "LIVER", "LUNG",
           "OESOPHAGUS", "OVARY", "PANCREAS", "PLACENTA", "PLEURA", "PROSTATE", "SALIVARY_GLAND",
           "SKIN", "SMALL_INTESTINE", "SOFT_TISSUE", "STOMACH", "THYROID", "UPPER_AERODIGESTIVE_TRACT",
           "URINARY_TRACT")


def get_local_thresholds(tissue):
    
    return None


def discretize_gene_expression(df_sample_rpmk, df_local_thresholds,
                               lb_global, ub_global,
                               col_name='rpkm', gene_list=None,
                               as_frame=True):

    confidences = [-1, 1, 2, 3]
    if gene_list is None:
        gene_list = df_sample_rpmk.index

    confidences_dict = {}
    for gene in gene_list:
        if gene not in df_local_thresholds.index:
            confidences_dict[gene] = 0
            continue

        rpkm = df_sample_rpmk[col_name][gene]
        local_threshold = df_local_thresholds[col_name][gene]
        thresholds = sorted([local_threshold, lb_global, ub_global, np.inf])
        for th, conf in zip(thresholds, confidences):
            if rpkm <= th:
                confidences_dict[gene] = conf
                break

    if as_frame:
        df_confidences = pd.DataFrame(data=list(confidences_dict.items()), columns=['hgnc_id', 'confidence'])
        df_confidences.set_index('hgnc_id', inplace=True)
        result = df_confidences
    else:
        result = confidences_dict

    return result
