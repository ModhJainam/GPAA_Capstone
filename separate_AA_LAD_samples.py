import pandas as pd

gpaa = pd.read_csv("GPAA_all_FPKM.csv")
meta = pd.read_csv("/project/gpaa/machine_learning/jainam_capstone/GPAA_samples_batches0-15_metadata_detail_grading_2022_02_22.csv")
LAD_meta = meta.loc[meta['tissue'].isin(['LAD44','LADLM1'])]
AA_meta = meta.loc[meta['tissue'].isin(['AA','AA18L', 'AA16L'])]
gpaa.replace(regex = ['L44', '-A16', '-A18'], value = ['LAD44', '-AA16', '-AA18'], inplace = True)
LAD = gpaa.loc[gpaa['sample_id'].isin(LAD_meta['UVA_sample_id'])]
AA = gpaa.loc[gpaa['sample_id'].isin(AA_meta['UVA_sample_id'])]
LAD.to_csv("GPAA_samples_LAD_batches0-15_FPKM.csv", index = False)
AA.to_csv("GPAA_samples_AA_batches0-15_FPKM.csv", index = False)
