import os
import glob
import shutil

fileList = [
    "MB_pPb8160_etaDijet_data_80*", "MB_pPb8160_etaDijet_data_90*", "MB_pPb8160_etaDijet_data_1*", "MB_pPb8160_etaDijet_data_2*", "MB_pPb8160_etaDijet_data_3*", "MB_pPb8160_etaDijet_data_400*", "MB_pPb8160_etaDijet_data_500*",
    "MB_pPb8160_etaDijet_pileupSyst_80*", "MB_pPb8160_etaDijet_pileupSyst_90*", "MB_pPb8160_etaDijet_pileupSyst_1*", "MB_pPb8160_etaDijet_pileupSyst_2*", "MB_pPb8160_etaDijet_pileupSyst_3*", "MB_pPb8160_etaDijet_pileupSyst_400*", "MB_pPb8160_etaDijet_pileupSyst_500*",
    "Jet60_pPb8160_etaDijet_*_100_110*", "Jet60_pPb8160_etaDijet_*_11*", "Jet60_pPb8160_etaDijet_*_12*", "Jet60_pPb8160_etaDijet_*_13*", "Jet60_pPb8160_etaDijet_*_14*", "Jet60_pPb8160_etaDijet_*_15*", "Jet60_pPb8160_etaDijet_*_16*", "Jet60_pPb8160_etaDijet_*_17*", "Jet60_pPb8160_etaDijet_*_18*", "Jet60_pPb8160_etaDijet_*_2*", "Jet60_pPb8160_etaDijet_*_3*", "Jet60_pPb8160_etaDijet_*_4*", "Jet60_pPb8160_etaDijet_*_5*", "Jet60_pPb8160_etaDijet_*_6*", "Jet60_pPb8160_etaDijet_*_7*",
    "Jet80_pPb8160_etaDijet_*_4*", "Jet80_pPb8160_etaDijet_*_5*", "Jet80_pPb8160_etaDijet_*_6*", "Jet80_pPb8160_etaDijet_*_7*", "Jet80_pPb8160_etaDijet_*_8*", "Jet80_pPb8160_etaDijet_*_9*",
    "Jet100_pPb8160_etaDijet_*_40_*", "Jet100_pPb8160_etaDijet_*_50_*", "Jet100_pPb8160_etaDijet_*_60_*", "Jet100_pPb8160_etaDijet_*_70*", "Jet100_pPb8160_etaDijet_*_80*", "Jet100_pPb8160_etaDijet_*_90*",
    "Jet100_pPb8160_etaDijet_*_100_*", "Jet100_pPb8160_etaDijet_*_110_*"
]

for pattern in fileList:
    for file in glob.glob(pattern):
        if (os.path.exists(file)):
            os.remove(file)
            print(f"{file} deleted")
        else:
            print(f"{file} does not exist")


# Copy pileup values from 150-1000 GeV pTave bin to the high-pt bins
pileup_original = "Jet100_pPb8160_etaDijet_pileupSyst_150_1000_lab.txt"
copy_pileup_syst = [ "Jet100_pPb8160_etaDijet_pileupSyst_150_160_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_160_180_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_180_200_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_200_240_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_240_300_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_300_400_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_400_500_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_500_1000_lab.txt", "Jet100_pPb8160_etaDijet_pileupSyst_300_1000_lab.txt",]
for f in copy_pileup_syst:
    shutil.copyfile(pileup_original, f)
    print(f"cp {pileup_original} {f}")