# Cell-Type Constellations in the Single-Cell Universe

Code accompanying **Cell-Type Constellations in the Single-Cell Universe**.
This repository contains one final source-processing entry point for every
PBMC study and the exact analysis used for Figures 2 and 3. Historical and
superseded scripts are intentionally excluded.

## Contents

- `scripts/sources/`: one final Python entry point per source study. Each entry
  downloads reviewed public inputs, applies SCimilarity v1.1 annotation, and
  calculates sample-level cell-type stars.
- `configs/`: exact final dataset-level processing configurations.
- `src/`: shared download, annotation, normalization, and star
  calculation implementation.
- `scripts/`: filtering, assembly, reference construction, alignment,
  and Figure 2/3 scripts.
- `metadata/`: source inventory and the ranked 2,000-gene panel.

Raw counts, the SCimilarity model, and generated matrices are not redistributed.
Source access terms remain those of the original repositories.

## Environment

Python 3.12 and a CUDA-capable PyTorch installation are recommended.

```bash
conda env create -f environment.yml
conda activate cell-type-constellations
```

Download the SCimilarity v1.1 model separately and place it at
`models/model_v1.1/`.

## Reproducing the star maps

Run a study entry point from the repository root, for example:

```bash
python scripts/sources/GSE149689.py
```

Run all entries in the table below to rebuild the full collection. Each script
is restart-safe: validated outputs are checksum-verified and reused. Downloaded
count carriers are removed only after derived outputs pass validation.

## Reproducing Figures 2 and 3

After source processing has completed:

```bash
python scripts/freeze_pbmc_66_cohort.py
python scripts/select_pbmc_66_within_map_variable_genes.py
python scripts/materialize_pbmc_66_star_matrix.py
python scripts/build_pbmc_66_reference_and_alignment.py
python scripts/build_pbmc_66_figure2.py
python scripts/build_pbmc_66_figure3.py
```

Expected final dimensions are 66 cell types, 342,559 stars in 14,230 maps,
17,900 genes present in at least 90% of stars, and 2,000 selected genes.

## Runtime and compute requirements

| Stage | Practical requirement | Time |
|---|---|---|
| All downloads, annotation, and star calculation | 1 CUDA GPU recommended; at least 64 GB RAM; up to ~1.5 TB temporary storage if many sources are staged | Several days to weeks for ~93.8 million cells; network and source format dominate |
| Cohort filtering and gene selection | 32–64 GB RAM; CPU | About 20 minutes |
| Raw and aligned matrix construction | 32–64 GB RAM; ~15 GB free storage | About 20–40 minutes |
| Figure 2 | CUDA GPU recommended; 32 GB RAM | About 2 minutes after matrices exist |
| Figure 3 | 8 CPU threads; 32 GB RAM | About 2 minutes after matrices exist |

The two principal float32 matrices are approximately 2.74 GB each. Exact timing
depends on network throughput, GPU model, filesystem speed, and source format.

## Source studies

The table lists every study, its constituent dataset identifiers, retained cell
observations, and the single final processing script released for that study.

| Study/source | Repository category | Dataset identifier(s) | Retained cells | Final script |
|---|---|---|---:|---|
| 03f821b4-87be-4ff4-b65a-b5fc00061da7 | CELLxGENE | 2a498ace-872a-4935-984b-1afa70fd9886 | 339,558 | `scripts/sources/03f821b4-87be-4ff4-b65a-b5fc00061da7.py` |
| 0a839c4b-10d0-4d64-9272-684c49a2c8ba | CELLxGENE, GEO-associated | 9dbab10c-118d-496b-966a-67f1763a6b7d<br>GSE158055_original_additional_views | 1,053,579 | `scripts/sources/0a839c4b-10d0-4d64-9272-684c49a2c8ba.py` |
| 29f92179-ca10-4309-a32b-d383d80347c1 | CELLxGENE | 01c93cf6-b695-4e30-a26e-121ae8b16a9e<br>50eb1e23-b8d4-4f76-a184-44e5541fa05a<br>79ef1959-a6b4-4cac-82ca-30feaec48df1<br>84230ea4-998d-4aa8-8456-81dd54ce23af<br>881fe679-c6e0-45a3-9427-c4e81be6921f<br>8a554710-08bc-4005-87cd-da9675bdc2e7<br>c874f155-9bf9-4928-b821-f52c876b3e48<br>cab0bc48-744c-461b-aaf4-7bf2cb7af00d<br>db59611b-42de-4035-93aa-1ed39f38b467<br>ea786a06-5855-48b7-80d7-0313a21a2044<br>ed9e9f96-4f08-49d2-bef5-b2c29adf3edc | 60,822 | `scripts/sources/29f92179-ca10-4309-a32b-d383d80347c1.py` |
| 3a518635-c062-4c66-a9bd-c4581a5e3c63 | CELLxGENE | 789ad837-9bf5-436c-9b60-9f1153546e4c | 2,029,835 | `scripts/sources/3a518635-c062-4c66-a9bd-c4581a5e3c63.py` |
| 40272c3b-4697-4bd4-ba3f-82fa96b9bf71 | Human Cell Atlas | HCA_eqtl45 | 25,183 | `scripts/sources/40272c3b-4697-4bd4-ba3f-82fa96b9bf71.py` |
| 436154da-bcf1-4130-9c8b-120ff9a888f2 | CELLxGENE | 218acb0f-9f2f-4f76-b90b-15a4b7c7f629 | 1,258,161 | `scripts/sources/436154da-bcf1-4130-9c8b-120ff9a888f2.py` |
| 4f99649a-2021-47a3-9908-2b83c546cbc9 | CELLxGENE | e365d9c0-f30c-48b5-9839-f93ee12f9013 | 404,171 | `scripts/sources/4f99649a-2021-47a3-9908-2b83c546cbc9.py` |
| 5542eeb0-96ef-4ab9-95ea-eb6abc178461 | CELLxGENE | 4a5b00e0-1ba3-4fd4-af89-d3512eb20720<br>d9f3a82c-24bb-4277-b715-8ec92a33dd32 | 1,011,797 | `scripts/sources/5542eeb0-96ef-4ab9-95ea-eb6abc178461.py` |
| 60a2676d-9f37-46cc-9b02-c7370a53be9c | CELLxGENE | ff4235bb-1cdf-4dcf-a13a-dc4f32244400 | 3,758,514 | `scripts/sources/60a2676d-9f37-46cc-9b02-c7370a53be9c.py` |
| 77f9d7e9-5675-49c3-abed-ce02f39eef1b | Author-hosted atlas resources, CELLxGENE | 0a8e3443-c3e2-4918-84db-0495657d9175<br>0c0c49a2-b96b-432f-98a4-1decf9929149<br>3cb646ab-d411-4e5b-89eb-082a786c953b<br>420d0f5d-b7ba-4b1e-88d0-44e9f733febd<br>46104f0b-9af5-466a-ae0f-56b8dc1969a2<br>47e00f43-6176-44a5-b4d1-f88a1920fcab<br>5fde5c9c-5b1e-4df5-982a-3f8e7635161f<br>AIFI_human_immune_health_atlas_full<br>aa32098f-f7b0-43f2-a00a-601d190aae94 | 14,071,262 | `scripts/sources/77f9d7e9-5675-49c3-abed-ce02f39eef1b.py` |
| 7d7cabfd-1d1f-40af-96b7-26a0825a306d | CELLxGENE | 01ad3cd7-3929-4654-84c0-6db05bd5fd59 | 462,362 | `scripts/sources/7d7cabfd-1d1f-40af-96b7-26a0825a306d.py` |
| 7f7fdf50-aa0e-4a6f-979c-53655cac9d5d | CELLxGENE | 0e5c75da-7bc9-4a2c-83fc-ef1e378bd5f9 | 481,007 | `scripts/sources/7f7fdf50-aa0e-4a6f-979c-53655cac9d5d.py` |
| 896f377c-8e88-463e-82b0-b2a5409d6fe4 | Author-hosted atlas resources | ABF300_immune_aging2026 | 987,997 | `scripts/sources/896f377c-8e88-463e-82b0-b2a5409d6fe4.py` |
| 8f126edf-5405-4731-8374-b5ce11f53e82 | CELLxGENE | ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3 | 834,709 | `scripts/sources/8f126edf-5405-4731-8374-b5ce11f53e82.py` |
| CIMA | Author-hosted atlas resources | CIMA_native_counts | 6,484,974 | `scripts/sources/CIMA.py` |
| E-MTAB-15719 | ArrayExpress/BioStudies | E-MTAB-15719_LC1-3<br>E-MTAB-15719_LC1_1-2<br>E-MTAB-15719_LC2<br>E-MTAB-15719_LC3 | 62,287 | `scripts/sources/E-MTAB-15719.py` |
| EGAS00001004571 | Human Cell Atlas | EGAS00001004571_PBMCs | 44,458 | `scripts/sources/EGAS00001004571.py` |
| GSE106540 | GEO-associated | GSE106540 | 2,186 | `scripts/sources/GSE106540.py` |
| GSE111014 | GEO-associated | GSE111014_genes001 | 40,294 | `scripts/sources/GSE111014.py` |
| GSE115189 | GEO-associated | GSE115189_genes001 | 3,364 | `scripts/sources/GSE115189.py` |
| GSE121267 | GEO-associated | GSE121267_genes001 | 3,148 | `scripts/sources/GSE121267.py` |
| GSE127471 | GEO-associated | GSE127471_genes001 | 1,801 | `scripts/sources/GSE127471.py` |
| GSE128243 | GEO-associated | GSE128243_genes001 | 13,790 | `scripts/sources/GSE128243.py` |
| GSE130117 | GEO-associated | GSE130117_genes001 | 2,863 | `scripts/sources/GSE130117.py` |
| GSE130157 | GEO-associated | GSE130157_genes001 | 69,598 | `scripts/sources/GSE130157.py` |
| GSE132065 | GEO-associated | GSE132065_clinical_RT<br>GSE132065_healthy_JULIA03 | 10,981 | `scripts/sources/GSE132065.py` |
| GSE132338 | GEO-associated | GSE132338_PBMC_component | 52,055 | `scripts/sources/GSE132338.py` |
| GSE132802 | GEO-associated | GSE132802_genes001<br>GSE132802_genes002 | 22,754 | `scripts/sources/GSE132802.py` |
| GSE134809 | GEO-associated | GSE134809_genes001 | 39,245 | `scripts/sources/GSE134809.py` |
| GSE135325 | GEO-associated | GSE135325_genes001 | 2,116 | `scripts/sources/GSE135325.py` |
| GSE135779 | GEO-associated | GSE135779 | 358,235 | `scripts/sources/GSE135779.py` |
| GSE136103 | GEO-associated | GSE136103 | 30,734 | `scripts/sources/GSE136103.py` |
| GSE138266 | GEO-associated | GSE138266 | 41,128 | `scripts/sources/GSE138266.py` |
| GSE139324 | GEO-associated | GSE139324 | 59,121 | `scripts/sources/GSE139324.py` |
| GSE140881 | GEO-associated | GSE140881_genes001<br>GSE140881_genes002<br>GSE140881_genes003 | 77,076 | `scripts/sources/GSE140881.py` |
| GSE142016 | GEO-associated | GSE142016_genes001 | 24,009 | `scripts/sources/GSE142016.py` |
| GSE145307 | GEO-associated | GSE145307_genes001 | 34,903 | `scripts/sources/GSE145307.py` |
| GSE145809 | GEO-associated | GSE145809_genes001 | 17,066 | `scripts/sources/GSE145809.py` |
| GSE147794 | GEO-associated | GSE147794_genes001 | 5,152 | `scripts/sources/GSE147794.py` |
| GSE149689 | GEO-associated | GSE149689 | 60,688 | `scripts/sources/GSE149689.py` |
| GSE150861 | GEO-associated | GSE150861 | 16,336 | `scripts/sources/GSE150861.py` |
| GSE151531 | GEO-associated | GSE151531 | 623 | `scripts/sources/GSE151531.py` |
| GSE153098 | GEO-associated | GSE153098_genes001 | 12,079 | `scripts/sources/GSE153098.py` |
| GSE154386 | GEO-associated | GSE154386 | 182,192 | `scripts/sources/GSE154386.py` |
| GSE155223 | GEO-associated | GSE155223_genes001 | 76,464 | `scripts/sources/GSE155223.py` |
| GSE155673 | GEO-associated | GSE155673_genes001 | 56,591 | `scripts/sources/GSE155673.py` |
| GSE155698 | GEO-associated | GSE155698 | 77,875 | `scripts/sources/GSE155698.py` |
| GSE156989 | GEO-associated | GSE156989_baseline_genes001 | 86,968 | `scripts/sources/GSE156989.py` |
| GSE157007 | GEO-associated | GSE157007 | 112,857 | `scripts/sources/GSE157007.py` |
| GSE157278 | GEO-associated | GSE157278 | 58,589 | `scripts/sources/GSE157278.py` |
| GSE157447 | GEO-associated | GSE157447 | 764 | `scripts/sources/GSE157447.py` |
| GSE157829 | GEO-associated | GSE157829_genes001<br>GSE157829_genes002 | 34,412 | `scripts/sources/GSE157829.py` |
| GSE158769 | GEO-associated | GSE158769 | 500,089 | `scripts/sources/GSE158769.py` |
| GSE160097 | GEO-associated | GSE160097_genes001 | 47,489 | `scripts/sources/GSE160097.py` |
| GSE161918 | GEO-associated | GSE161918_native_reassembled | 343,798 | `scripts/sources/GSE161918.py` |
| GSE162117 | GEO-associated | GSE162117_source_donor_days | 113,942 | `scripts/sources/GSE162117.py` |
| GSE162498 | GEO-associated | GSE162498_genes001<br>GSE162498_genes002<br>GSE162498_genes003 | 28,790 | `scripts/sources/GSE162498.py` |
| GSE162577 | GEO-associated | GSE162577 | 12,207 | `scripts/sources/GSE162577.py` |
| GSE162632 | GEO-associated | GSE162632_genes001 | 124,917 | `scripts/sources/GSE162632.py` |
| GSE162806 | GEO-associated | GSE162806_genes001<br>GSE162806_genes002<br>GSE162806_genes003<br>GSE162806_genes004<br>GSE162806_genes005<br>GSE162806_genes006<br>GSE162806_genes007<br>GSE162806_genes008<br>GSE162806_genes009<br>GSE162806_genes010 | 176,881 | `scripts/sources/GSE162806.py` |
| GSE163121 | GEO-associated | GSE163121_genes001 | 25,001 | `scripts/sources/GSE163121.py` |
| GSE163314 | GEO-associated | GSE163314_genes001 | 50,515 | `scripts/sources/GSE163314.py` |
| GSE164378 | GEO-associated | GSE164378_3P | 161,764 | `scripts/sources/GSE164378.py` |
| GSE164379 | GEO-associated | GSE164379_genes001 | 31,658 | `scripts/sources/GSE164379.py` |
| GSE165080 | GEO-associated | GSE165080 | 200,528 | `scripts/sources/GSE165080.py` |
| GSE165087 | GEO-associated | GSE165087_genes001 | 383,638 | `scripts/sources/GSE165087.py` |
| GSE165816 | GEO-associated | GSE165816_G19<br>GSE165816_G20<br>GSE165816_G21<br>GSE165816_G22<br>GSE165816_G25<br>GSE165816_G26<br>GSE165816_G27<br>GSE165816_G30<br>GSE165816_G37<br>GSE165816_G47 | 43,288 | `scripts/sources/GSE165816.py` |
| GSE166181 | GEO-associated | GSE166181_genes001 | 35,108 | `scripts/sources/GSE166181.py` |
| GSE166352 | GEO-associated | GSE166352_genes001 | 38,762 | `scripts/sources/GSE166352.py` |
| GSE166489 | GEO-associated | GSE166489_genes001 | 162,871 | `scripts/sources/GSE166489.py` |
| GSE166992 | GEO-associated | GSE166992_genes001 | 60,077 | `scripts/sources/GSE166992.py` |
| GSE167029 | GEO-associated | GSE167029_genes001 | 198,861 | `scripts/sources/GSE167029.py` |
| GSE167159 | GEO-associated | GSE167159_genes001 | 42,843 | `scripts/sources/GSE167159.py` |
| GSE167363 | GEO-associated | GSE167363 | 59,754 | `scripts/sources/GSE167363.py` |
| GSE167976 | GEO-associated | GSE167976 | 55,656 | `scripts/sources/GSE167976.py` |
| GSE168732 | GEO-associated | GSE168732 | 103,362 | `scripts/sources/GSE168732.py` |
| GSE169246 | GEO-associated | GSE169246_gene_rna | 322,246 | `scripts/sources/GSE169246.py` |
| GSE171555 | GEO-associated | GSE171555 | 513,561 | `scripts/sources/GSE171555.py` |
| GSE171811 | GEO-associated | GSE171811_38858e461760<br>GSE171811_a7bfe9ffe578 | 30,293 | `scripts/sources/GSE171811.py` |
| GSE172179 | GEO-associated | GSE172179_genes001 | 7,107 | `scripts/sources/GSE172179.py` |
| GSE173303 | GEO-associated | GSE173303 | 59,254 | `scripts/sources/GSE173303.py` |
| GSE173932 | GEO-associated | GSE173932_genes001 | 16,129 | `scripts/sources/GSE173932.py` |
| GSE174461 | GEO-associated | GSE174461_genes001 | 5,838 | `scripts/sources/GSE174461.py` |
| GSE175450 | GEO-associated | GSE175450 | 58,823 | `scripts/sources/GSE175450.py` |
| GSE175930 | GEO-associated | GSE175930_genes001 | 30,352 | `scripts/sources/GSE175930.py` |
| GSE178429 | GEO-associated | GSE178429 | 22,436 | `scripts/sources/GSE178429.py` |
| GSE178882 | GEO-associated | GSE178882_genes001 | 109,650 | `scripts/sources/GSE178882.py` |
| GSE180045 | GEO-associated | GSE180045_genes001<br>GSE180045_genes002 | 101,116 | `scripts/sources/GSE180045.py` |
| GSE180578 | GEO-associated | GSE180578_genes001<br>GSE180578_genes002 | 85,901 | `scripts/sources/GSE180578.py` |
| GSE181279 | GEO-associated | GSE181279 | 36,751 | `scripts/sources/GSE181279.py` |
| GSE181543 | GEO-associated | GSE181543_genes001 | 12,547 | `scripts/sources/GSE181543.py` |
| GSE181897 | GEO-associated | GSE181897_nativeRNA | 135,676 | `scripts/sources/GSE181897.py` |
| GSE182123 | GEO-associated | GSE182123_genes001 | 14,115 | `scripts/sources/GSE182123.py` |
| GSE182536 | GEO-associated | GSE182536_genes001 | 21,524 | `scripts/sources/GSE182536.py` |
| GSE183953 | GEO-associated | GSE183953_genes001 | 59,519 | `scripts/sources/GSE183953.py` |
| GSE184031 | GEO-associated | GSE184031_genes001 | 5,849 | `scripts/sources/GSE184031.py` |
| GSE184241 | GEO-associated | GSE184241 | 1,710 | `scripts/sources/GSE184241.py` |
| GSE185712 | GEO-associated | GSE185712 | 14,010 | `scripts/sources/GSE185712.py` |
| GSE185857 | GEO-associated | GSE185857_genes001 | 60,222 | `scripts/sources/GSE185857.py` |
| GSE188172 | GEO-associated | GSE188172_GSM5671235<br>GSE188172_GSM5671236<br>GSE188172_GSM5671237<br>GSE188172_GSM5671238<br>GSE188172_GSM5671239<br>GSE188172_GSM5671240<br>GSE188172_GSM5671241<br>GSE188172_GSM5671242<br>GSE188172_GSM5671243<br>GSE188172_GSM5671244<br>GSE188172_GSM5671245<br>GSE188172_GSM5671246<br>GSE188172_GSM5671247<br>GSE188172_GSM5671248<br>GSE188172_GSM5671249<br>GSE188172_GSM5671250<br>GSE188172_GSM5671251<br>GSE188172_GSM5671252<br>GSE188172_GSM5671253<br>GSE188172_GSM5671254<br>GSE188172_GSM5671255<br>GSE188172_GSM5671256<br>GSE188172_GSM5671257<br>GSE188172_GSM5671258<br>GSE188172_GSM5671259<br>GSE188172_GSM5671260<br>GSE188172_GSM5671261<br>GSE188172_GSM5671262<br>GSE188172_GSM5671263<br>GSE188172_GSM5671264<br>GSE188172_GSM5671265<br>GSE188172_GSM5671266<br>GSE188172_GSM5671267<br>GSE188172_GSM5671268<br>GSE188172_GSM5671269<br>GSE188172_GSM5671270<br>GSE188172_GSM5671271<br>GSE188172_GSM5671272<br>GSE188172_GSM5671273<br>GSE188172_GSM5671274<br>GSE188172_GSM5671275<br>GSE188172_GSM5671276<br>GSE188172_GSM5671277<br>GSE188172_GSM5671278<br>GSE188172_GSM5671279<br>GSE188172_GSM5671280<br>GSE188172_GSM5671281<br>GSE188172_GSM5671282<br>GSE188172_GSM5671283<br>GSE188172_GSM5671284<br>GSE188172_GSM5671285<br>GSE188172_GSM5671286<br>GSE188172_GSM5671287<br>GSE188172_GSM5671288<br>GSE188172_GSM5671289<br>GSE188172_GSM5671290<br>GSE188172_GSM5671291<br>GSE188172_GSM5671292<br>GSE188172_GSM5671293<br>GSE188172_GSM5671294<br>GSE188172_GSM5671295<br>GSE188172_GSM5671296<br>GSE188172_GSM5671297<br>GSE188172_GSM5671298<br>GSE188172_GSM5671299<br>GSE188172_GSM5671300<br>GSE188172_GSM5671301<br>GSE188172_GSM5671302<br>GSE188172_GSM5671303<br>GSE188172_GSM5671304<br>GSE188172_GSM5671305<br>GSE188172_GSM5671306<br>GSE188172_GSM5671307<br>GSE188172_GSM5671308<br>GSE188172_GSM5671309<br>GSE188172_GSM5671310<br>GSE188172_GSM5671311<br>GSE188172_GSM5671312<br>GSE188172_GSM5671313<br>GSE188172_GSM5671314<br>GSE188172_GSM5671315<br>GSE188172_GSM5671316<br>GSE188172_GSM5671317<br>GSE188172_GSM5671318<br>GSE188172_GSM5671319<br>GSE188172_GSM5671320<br>GSE188172_GSM5671321<br>GSE188172_GSM5671322<br>GSE188172_GSM5671323<br>GSE188172_GSM5671324<br>GSE188172_GSM5671325<br>GSE188172_GSM5671326<br>GSE188172_GSM5671327<br>GSE188172_GSM5671328<br>GSE188172_GSM5671329<br>GSE188172_GSM5671330<br>GSE188172_GSM5671331<br>GSE188172_GSM5671332<br>GSE188172_GSM5671333<br>GSE188172_GSM5671334<br>GSE188172_GSM5671335<br>GSE188172_GSM5671336<br>GSE188172_GSM5671337<br>GSE188172_GSM5671338<br>GSE188172_GSM5671339<br>GSE188172_GSM5671340<br>GSE188172_GSM5671341<br>GSE188172_GSM5671342<br>GSE188172_GSM5671343<br>GSE188172_GSM5671344<br>GSE188172_GSM5671345<br>GSE188172_GSM5671346<br>GSE188172_GSM5671347<br>GSE188172_GSM5671348<br>GSE188172_GSM5671349<br>GSE188172_GSM5671350<br>GSE188172_GSM5671351<br>GSE188172_GSM5671352<br>GSE188172_GSM5671353<br>GSE188172_GSM5671354<br>GSE188172_GSM5671355<br>GSE188172_GSM5671356<br>GSE188172_GSM5671357<br>GSE188172_GSM5671358<br>GSE188172_GSM5671359<br>GSE188172_GSM5671360<br>GSE188172_GSM5671361<br>GSE188172_GSM5671362<br>GSE188172_GSM5671363<br>GSE188172_GSM5671364<br>GSE188172_GSM5671365<br>GSE188172_GSM5671366<br>GSE188172_GSM5671367<br>GSE188172_GSM5671368<br>GSE188172_GSM5671369<br>GSE188172_GSM5671370<br>GSE188172_GSM5671371<br>GSE188172_GSM5671372<br>GSE188172_GSM5671373<br>GSE188172_GSM5671374<br>GSE188172_GSM5671375<br>GSE188172_GSM5671376<br>GSE188172_GSM5671377<br>GSE188172_GSM5671378<br>GSE188172_GSM5671379<br>GSE188172_GSM5671383<br>GSE188172_GSM5671384<br>GSE188172_GSM5671385<br>GSE188172_GSM5671386<br>GSE188172_GSM5671387<br>GSE188172_GSM5671388<br>GSE188172_GSM5671389<br>GSE188172_GSM5671390<br>GSE188172_GSM5671391<br>GSE188172_GSM5671392<br>GSE188172_GSM5671393<br>GSE188172_GSM5671394<br>GSE188172_GSM5671395<br>GSE188172_GSM5671396<br>GSE188172_GSM5671397<br>GSE188172_GSM5671398<br>GSE188172_GSM5671399<br>GSE188172_GSM5671400<br>GSE188172_GSM5671401<br>GSE188172_GSM5671402<br>GSE188172_GSM5671403<br>GSE188172_GSM5671404<br>GSE188172_GSM5671405 | 244,091 | `scripts/sources/GSE188172.py` |
| GSE188620 | GEO-associated | GSE188620 | 151,740 | `scripts/sources/GSE188620.py` |
| GSE189050 | GEO-associated | GSE189050 | 22,437 | `scripts/sources/GSE189050.py` |
| GSE189125 | GEO-associated | GSE189125 | 27,131 | `scripts/sources/GSE189125.py` |
| GSE189188 | GEO-associated | GSE189188_genes001 | 104,576 | `scripts/sources/GSE189188.py` |
| GSE190510 | GEO-associated | GSE190510 | 58,425 | `scripts/sources/GSE190510.py` |
| GSE190760 | GEO-associated | GSE190760_genes001 | 280,281 | `scripts/sources/GSE190760.py` |
| GSE190992 | GEO-associated | GSE190992_native_ptid2<br>GSE190992_native_ptid4<br>GSE190992_native_ptid5<br>GSE190992_native_ptid6 | 498,399 | `scripts/sources/GSE190992.py` |
| GSE191332 | GEO-associated | GSE191332_visit1<br>GSE191332_visit2<br>GSE191332_visit3 | 13,034 | `scripts/sources/GSE191332.py` |
| GSE192391 | GEO-associated | GSE192391 | 216,591 | `scripts/sources/GSE192391.py` |
| GSE193220 | GEO-associated | GSE193220_genes001 | 24,711 | `scripts/sources/GSE193220.py` |
| GSE193525 | GEO-associated | GSE193525 | 171,241 | `scripts/sources/GSE193525.py` |
| GSE193770 | GEO-associated | GSE193770_genes001 | 67,295 | `scripts/sources/GSE193770.py` |
| GSE194078 | GEO-associated | GSE194078_genes001 | 109,068 | `scripts/sources/GSE194078.py` |
| GSE194245 | GEO-associated | GSE194245_genes001 | 10,463 | `scripts/sources/GSE194245.py` |
| GSE194315 | GEO-associated | GSE194315 | 271,883 | `scripts/sources/GSE194315.py` |
| GSE197268 | GEO-associated | GSE197268_genes001 | 208,236 | `scripts/sources/GSE197268.py` |
| GSE197658 | GEO-associated | GSE197658_RNA_18A<br>GSE197658_RNA_18D<br>GSE197658_RNA_18N<br>GSE197658_RNA_23A<br>GSE197658_RNA_23D<br>GSE197658_RNA_23N<br>GSE197658_RNA_30A<br>GSE197658_RNA_30D<br>GSE197658_RNA_30N<br>GSE197658_RNA_32A<br>GSE197658_RNA_32D<br>GSE197658_RNA_32N<br>GSE197658_RNA_35A<br>GSE197658_RNA_35D<br>GSE197658_RNA_35N<br>GSE197658_RNA_36A<br>GSE197658_RNA_36D<br>GSE197658_RNA_36N<br>GSE197658_RNA_37A<br>GSE197658_RNA_37D<br>GSE197658_RNA_37N<br>GSE197658_RNA_38A<br>GSE197658_RNA_38D<br>GSE197658_RNA_38N<br>GSE197658_RNA_39A<br>GSE197658_RNA_39D<br>GSE197658_RNA_39N<br>GSE197658_RNA_40A<br>GSE197658_RNA_40D<br>GSE197658_RNA_40N | 203,170 | `scripts/sources/GSE197658.py` |
| GSE198339 | GEO-associated | GSE198339 | 9,366 | `scripts/sources/GSE198339.py` |
| GSE198616 | GEO-associated | GSE198616_genes001<br>GSE198616_genes002<br>GSE198616_monocyte_genes001 | 88,976 | `scripts/sources/GSE198616.py` |
| GSE198891 | GEO-associated | GSE198891_genes001 | 31,768 | `scripts/sources/GSE198891.py` |
| GSE199445 | GEO-associated | GSE199445 | 75,248 | `scripts/sources/GSE199445.py` |
| GSE199994 | GEO-associated | GSE199994 | 96,954 | `scripts/sources/GSE199994.py` |
| GSE200107 | GEO-associated | GSE200107_genes001<br>GSE200107_genes002 | 126,095 | `scripts/sources/GSE200107.py` |
| GSE200996 | GEO-associated | GSE200996 | 145,584 | `scripts/sources/GSE200996.py` |
| GSE201425 | GEO-associated | GSE201425_genes001 | 5,752 | `scripts/sources/GSE201425.py` |
| GSE201534 | GEO-associated | GSE201534_genes001 | 73,224 | `scripts/sources/GSE201534.py` |
| GSE204716 | GEO-associated | GSE204716_genes001 | 92,948 | `scripts/sources/GSE204716.py` |
| GSE205704 | GEO-associated | GSE205704_genes001 | 26,144 | `scripts/sources/GSE205704.py` |
| GSE206265 | GEO-associated | GSE206265_participants | 609,492 | `scripts/sources/GSE206265.py` |
| GSE206283 | GEO-associated | GSE206283_genes001 | 458,709 | `scripts/sources/GSE206283.py` |
| GSE206295 | GEO-associated | GSE206295 | 32,870 | `scripts/sources/GSE206295.py` |
| GSE206298 | GEO-associated | GSE206298 | 118,817 | `scripts/sources/GSE206298.py` |
| GSE207192 | GEO-associated | GSE207192 | 55,138 | `scripts/sources/GSE207192.py` |
| GSE207475 | GEO-associated | GSE207475 | 81,012 | `scripts/sources/GSE207475.py` |
| GSE207633 | GEO-associated | GSE207633 | 199,220 | `scripts/sources/GSE207633.py` |
| GSE207679 | GEO-associated | GSE207679_genes001 | 135,998 | `scripts/sources/GSE207679.py` |
| GSE207935 | GEO-associated | GSE207935 | 140,575 | `scripts/sources/GSE207935.py` |
| GSE208337 | GEO-associated | GSE208337_genes001 | 99,616 | `scripts/sources/GSE208337.py` |
| GSE210065 | GEO-associated | GSE210065_PBMC | 244,403 | `scripts/sources/GSE210065.py` |
| GSE210433 | GEO-associated | GSE210433<br>GSE210433_multiome_genes001 | 123,834 | `scripts/sources/GSE210433.py` |
| GSE211068 | GEO-associated | GSE211068_genes001 | 18,041 | `scripts/sources/GSE211068.py` |
| GSE211504 | GEO-associated | GSE211504_genes001 | 38,071 | `scripts/sources/GSE211504.py` |
| GSE211560 | GEO-associated | GSE211560_genes001 | 57,443 | `scripts/sources/GSE211560.py` |
| GSE212038 | GEO-associated | GSE212038_genes001 | 19,637 | `scripts/sources/GSE212038.py` |
| GSE212092 | GEO-associated | GSE212092_genes001 | 76,002 | `scripts/sources/GSE212092.py` |
| GSE212217 | GEO-associated | GSE212217_genes001 | 298,683 | `scripts/sources/GSE212217.py` |
| GSE212479 | GEO-associated | GSE212479 | 5,625 | `scripts/sources/GSE212479.py` |
| GSE212584 | GEO-associated | GSE212584_genes001<br>GSE212584_genes002 | 115,590 | `scripts/sources/GSE212584.py` |
| GSE213095 | GEO-associated | GSE213095_genes001<br>GSE213095_genes002 | 91,347 | `scripts/sources/GSE213095.py` |
| GSE213243 | GEO-associated | GSE213243_genes001 | 20,220 | `scripts/sources/GSE213243.py` |
| GSE213516 | GEO-associated | GSE213516_genes001 | 136,673 | `scripts/sources/GSE213516.py` |
| GSE214283 | GEO-associated | GSE214283_genes001 | 503,900 | `scripts/sources/GSE214283.py` |
| GSE214546 | GEO-associated | GSE214546_pediatric_native | 140,776 | `scripts/sources/GSE214546.py` |
| GSE214865 | GEO-associated | GSE214865 | 109,720 | `scripts/sources/GSE214865.py` |
| GSE214985 | GEO-associated | GSE214985_genes001 | 113,001 | `scripts/sources/GSE214985.py` |
| GSE215219 | GEO-associated | GSE215219 | 66,948 | `scripts/sources/GSE215219.py` |
| GSE216329 | GEO-associated | GSE216329_genes001<br>GSE216329_genes002<br>GSE216329_genes003<br>GSE216329_genes004 | 164,370 | `scripts/sources/GSE216329.py` |
| GSE216486 | GEO-associated | GSE216486_genes001 | 6,952 | `scripts/sources/GSE216486.py` |
| GSE216548 | GEO-associated | GSE216548_genes001 | 73,949 | `scripts/sources/GSE216548.py` |
| GSE216762 | GEO-associated | GSE216762_genes001 | 14,398 | `scripts/sources/GSE216762.py` |
| GSE217561 | GEO-associated | GSE217561 | 124,630 | `scripts/sources/GSE217561.py` |
| GSE217906 | GEO-associated | GSE217906_genes001 | 46,018 | `scripts/sources/GSE217906.py` |
| GSE217930 | GEO-associated | GSE217930_genes001<br>GSE217930_genes002 | 109,539 | `scripts/sources/GSE217930.py` |
| GSE219086 | GEO-associated | GSE219086 | 56,526 | `scripts/sources/GSE219086.py` |
| GSE220699 | GEO-associated | GSE220699_genes001 | 119,694 | `scripts/sources/GSE220699.py` |
| GSE220764 | GEO-associated | GSE220764_genes001<br>GSE220764_genes002 | 114,430 | `scripts/sources/GSE220764.py` |
| GSE220790 | GEO-associated | GSE220790_genes001 | 65,206 | `scripts/sources/GSE220790.py` |
| GSE220969 | GEO-associated | GSE220969_genes001 | 157,013 | `scripts/sources/GSE220969.py` |
| GSE222200 | GEO-associated | GSE222200_genes001 | 5,852 | `scripts/sources/GSE222200.py` |
| GSE222647 | GEO-associated | GSE222647 | 177,582 | `scripts/sources/GSE222647.py` |
| GSE223138 | GEO-associated | GSE223138_genes001 | 63,671 | `scripts/sources/GSE223138.py` |
| GSE224198 | GEO-associated | GSE224198_genes001 | 306,380 | `scripts/sources/GSE224198.py` |
| GSE226039 | GEO-associated | GSE226039_genes001 | 13,186 | `scripts/sources/GSE226039.py` |
| GSE226099 | GEO-associated | GSE226099_genes001 | 196,283 | `scripts/sources/GSE226099.py` |
| GSE226328 | GEO-associated | GSE226328_genes001 | 16,628 | `scripts/sources/GSE226328.py` |
| GSE226557 | GEO-associated | GSE226557_genes001 | 19,289 | `scripts/sources/GSE226557.py` |
| GSE226598 | GEO-associated | GSE226598_genes001 | 111,954 | `scripts/sources/GSE226598.py` |
| GSE227001 | GEO-associated | GSE227001_genes001 | 131,755 | `scripts/sources/GSE227001.py` |
| GSE227835 | GEO-associated | GSE227835_genes001 | 427,810 | `scripts/sources/GSE227835.py` |
| GSE227954 | GEO-associated | GSE227954_genes001 | 17,517 | `scripts/sources/GSE227954.py` |
| GSE228279 | GEO-associated | GSE228279_genes001<br>GSE228279_genes002 | 18,722 | `scripts/sources/GSE228279.py` |
| GSE228597 | GEO-associated | GSE228597_controls_genes001 | 114,124 | `scripts/sources/GSE228597.py` |
| GSE228759 | GEO-associated | GSE228759_genes001 | 26,885 | `scripts/sources/GSE228759.py` |
| GSE229173 | GEO-associated | GSE229173_genes001<br>GSE229173_genes002<br>GSE229173_genes003<br>GSE229173_genes004<br>GSE229173_genes005<br>GSE229173_genes006 | 56,482 | `scripts/sources/GSE229173.py` |
| GSE229431 | GEO-associated | GSE229431_genes001 | 59,024 | `scripts/sources/GSE229431.py` |
| GSE229733 | GEO-associated | GSE229733 | 319,962 | `scripts/sources/GSE229733.py` |
| GSE230227 | GEO-associated | GSE230227_RNA | 182,456 | `scripts/sources/GSE230227.py` |
| GSE231794 | GEO-associated | GSE231794_genes001<br>GSE231794_genes002<br>GSE231794_genes003<br>GSE231794_genes004<br>GSE231794_genes005<br>GSE231794_genes006<br>GSE231794_genes007<br>GSE231794_genes008 | 133,501 | `scripts/sources/GSE231794.py` |
| GSE231906 | GEO-associated | GSE231906_genes001<br>GSE231906_genes002 | 246,784 | `scripts/sources/GSE231906.py` |
| GSE232186 | GEO-associated | GSE232186_genes001 | 35,513 | `scripts/sources/GSE232186.py` |
| GSE233046 | GEO-associated | GSE233046_genes001 | 110,155 | `scripts/sources/GSE233046.py` |
| GSE233275 | GEO-associated | GSE233275 | 79,590 | `scripts/sources/GSE233275.py` |
| GSE233276 | GEO-associated | GSE233276_genes001 | 66,649 | `scripts/sources/GSE233276.py` |
| GSE233304 | GEO-associated | GSE233304_genes001 | 35,315 | `scripts/sources/GSE233304.py` |
| GSE233321 | GEO-associated | GSE233321_469ab5415668<br>GSE233321_d8179ef6dba8 | 637,936 | `scripts/sources/GSE233321.py` |
| GSE233844 | GEO-associated | GSE233844 | 166,233 | `scripts/sources/GSE233844.py` |
| GSE233917 | GEO-associated | GSE233917_genes001 | 20,559 | `scripts/sources/GSE233917.py` |
| GSE234261 | GEO-associated | GSE234261_genes001<br>GSE234261_genes002 | 139,881 | `scripts/sources/GSE234261.py` |
| GSE234505 | GEO-associated | GSE234505 | 75,452 | `scripts/sources/GSE234505.py` |
| GSE234521 | GEO-associated | GSE234521 | 26,822 | `scripts/sources/GSE234521.py` |
| GSE235063 | GEO-associated | GSE235063_genes001 | 85,117 | `scripts/sources/GSE235063.py` |
| GSE235091 | GEO-associated | GSE235091_genes001 | 151,422 | `scripts/sources/GSE235091.py` |
| GSE235512 | GEO-associated | GSE235512_genes001 | 17,116 | `scripts/sources/GSE235512.py` |
| GSE235857 | GEO-associated | GSE235857 | 91,770 | `scripts/sources/GSE235857.py` |
| GSE236581 | GEO-associated | GSE236581 | 416,283 | `scripts/sources/GSE236581.py` |
| GSE237579 | GEO-associated | GSE237579_genes001 | 96,149 | `scripts/sources/GSE237579.py` |
| GSE237611 | GEO-associated | GSE237611_genes001 | 42,077 | `scripts/sources/GSE237611.py` |
| GSE238130 | GEO-associated | GSE238130_genes001 | 198,325 | `scripts/sources/GSE238130.py` |
| GSE239293 | GEO-associated | GSE239293_genes001 | 155 | `scripts/sources/GSE239293.py` |
| GSE239450 | GEO-associated | GSM7665278_OB7-CTRL_GEX<br>GSM7665279_705385-SEV_GEX<br>GSM7665280_803763-ASX_GEX<br>GSM7665281_324058-ASX_GEX<br>GSM7665282_644394-CTRL_GEX<br>GSM7665283_915348-SEV_GEX<br>GSM7665284_729106-CTRL_GEX<br>GSM7665285_105199-ASX_GEX<br>GSM7665286_101607-SEV_GEX | 120,448 | `scripts/sources/GSE239450.py` |
| GSE239626 | GEO-associated | GSE239626 | 40,830 | `scripts/sources/GSE239626.py` |
| GSE241096 | GEO-associated | GSE241096_genes001 | 63,730 | `scripts/sources/GSE241096.py` |
| GSE241492 | GEO-associated | GSE241492_genes001 | 53,867 | `scripts/sources/GSE241492.py` |
| GSE242780 | GEO-associated | GSE242780_genes001 | 45,961 | `scripts/sources/GSE242780.py` |
| GSE243002 | GEO-associated | GSE243002_genes001<br>GSE243002_genes002 | 76,476 | `scripts/sources/GSE243002.py` |
| GSE243061 | GEO-associated | GSE243061 | 43,320 | `scripts/sources/GSE243061.py` |
| GSE243572 | GEO-associated | GSE243572_genes001 | 111,305 | `scripts/sources/GSE243572.py` |
| GSE243629 | GEO-associated | GSE243629 | 102,826 | `scripts/sources/GSE243629.py` |
| GSE243905 | GEO-associated | GSE243905_genes001<br>GSE243905_genes002 | 237,127 | `scripts/sources/GSE243905.py` |
| GSE244173 | GEO-associated | GSE244173_genes001 | 99,249 | `scripts/sources/GSE244173.py` |
| GSE244263 | GEO-associated | GSE244263_genes001 | 340,603 | `scripts/sources/GSE244263.py` |
| GSE244515 | GEO-associated | GSE244515_genes001<br>GSE244515_genes002 | 283,358 | `scripts/sources/GSE244515.py` |
| GSE244589 | GEO-associated | GSE244589_genes001 | 83,021 | `scripts/sources/GSE244589.py` |
| GSE244793 | GEO-associated | GSE244793_genes001 | 166,332 | `scripts/sources/GSE244793.py` |
| GSE244831 | GEO-associated | GSE244831 | 73,905 | `scripts/sources/GSE244831.py` |
| GSE245131 | GEO-associated | GSE245131_genes001 | 50,083 | `scripts/sources/GSE245131.py` |
| GSE247277 | GEO-associated | GSE247277_genes001 | 61,143 | `scripts/sources/GSE247277.py` |
| GSE247322 | GEO-associated | GSE247322_genes001<br>GSE247322_genes002 | 172,146 | `scripts/sources/GSE247322.py` |
| GSE247824 | GEO-associated | GSE247824_unmixed_genes001 | 8,046 | `scripts/sources/GSE247824.py` |
| GSE247917 | GEO-associated | GSE247917_original_eligible | 180,934 | `scripts/sources/GSE247917.py` |
| GSE248256 | GEO-associated | GSE248256 | 106,924 | `scripts/sources/GSE248256.py` |
| GSE248284 | GEO-associated | GSE248284_genes001<br>GSE248284_genes002 | 47,921 | `scripts/sources/GSE248284.py` |
| GSE248556 | GEO-associated | GSE248556 | 152,236 | `scripts/sources/GSE248556.py` |
| GSE248685 | GEO-associated | GSE248685_genes001 | 38,220 | `scripts/sources/GSE248685.py` |
| GSE248689 | GEO-associated | GSE248689 | 37,416 | `scripts/sources/GSE248689.py` |
| GSE248693 | GEO-associated | GSE248693 | 144,309 | `scripts/sources/GSE248693.py` |
| GSE248695 | GEO-associated | GSE248695 | 206,314 | `scripts/sources/GSE248695.py` |
| GSE249131 | GEO-associated | GSE249131 | 184,132 | `scripts/sources/GSE249131.py` |
| GSE249313 | GEO-associated | GSE249313_genes001 | 320,865 | `scripts/sources/GSE249313.py` |
| GSE249513 | GEO-associated | GSE249513_genes001 | 69,868 | `scripts/sources/GSE249513.py` |
| GSE249584 | GEO-associated | GSE249584 | 55,557 | `scripts/sources/GSE249584.py` |
| GSE249898 | GEO-associated | GSE249898_genes001 | 306,666 | `scripts/sources/GSE249898.py` |
| GSE250024 | GEO-associated | GSE250024_genes001 | 30,920 | `scripts/sources/GSE250024.py` |
| GSE250488 | GEO-associated | GSE250488 | 16,544 | `scripts/sources/GSE250488.py` |
| GSE251696 | GEO-associated | GSE251696_genes001 | 56,567 | `scripts/sources/GSE251696.py` |
| GSE252273 | GEO-associated | GSE252273_verified_post | 32,931 | `scripts/sources/GSE252273.py` |
| GSE252572 | GEO-associated | GSE252572_genes001 | 51,848 | `scripts/sources/GSE252572.py` |
| GSE253173 | GEO-associated | GSE253173 | 385,038 | `scripts/sources/GSE253173.py` |
| GSE253279 | GEO-associated | GSE253279_genes001 | 207,363 | `scripts/sources/GSE253279.py` |
| GSE253485 | GEO-associated | GSE253485_genes001 | 11,254 | `scripts/sources/GSE253485.py` |
| GSE253568 | GEO-associated | GSE253568_genes001 | 105,775 | `scripts/sources/GSE253568.py` |
| GSE253633 | GEO-associated | GSE253633_genes001<br>GSE253633_genes002 | 50,730 | `scripts/sources/GSE253633.py` |
| GSE254176 | GEO-associated | GSE254176 | 118,541 | `scripts/sources/GSE254176.py` |
| GSE254249 | GEO-associated | GSE254249 | 372,346 | `scripts/sources/GSE254249.py` |
| GSE255566 | GEO-associated | GSE255566 | 38,386 | `scripts/sources/GSE255566.py` |
| GSE256089 | GEO-associated | GSE256089 | 140,171 | `scripts/sources/GSE256089.py` |
| GSE256291 | GEO-associated | GSE256291_genes001 | 24,538 | `scripts/sources/GSE256291.py` |
| GSE260763 | GEO-associated | GSE260763_genes001 | 259,907 | `scripts/sources/GSE260763.py` |
| GSE261140 | GEO-associated | GSE261140_genes001 | 16,622 | `scripts/sources/GSE261140.py` |
| GSE262146 | GEO-associated | GSE262146_genes001 | 125,734 | `scripts/sources/GSE262146.py` |
| GSE262240 | GEO-associated | GSE262240_genes001 | 88,588 | `scripts/sources/GSE262240.py` |
| GSE262861 | GEO-associated | GSE262861_GSM8180643<br>GSE262861_GSM8180644<br>GSE262861_GSM8180645<br>GSE262861_GSM8180646<br>GSE262861_GSM8180647<br>GSE262861_GSM8180648<br>GSE262861_GSM8180649<br>GSE262861_GSM8180650<br>GSE262861_GSM8180651<br>GSE262861_GSM8180652<br>GSE262861_GSM8180653 | 103,425 | `scripts/sources/GSE262861.py` |
| GSE263817 | GEO-associated | GSE263817_genes001 | 60,938 | `scripts/sources/GSE263817.py` |
| GSE263931 | GEO-associated | GSE263931_genes001 | 102,186 | `scripts/sources/GSE263931.py` |
| GSE263995 | GEO-associated | GSE263995_genes001 | 33,212 | `scripts/sources/GSE263995.py` |
| GSE264196 | GEO-associated | GSE264196_genes001 | 90,725 | `scripts/sources/GSE264196.py` |
| GSE264489 | GEO-associated | GSE264489 | 35,070 | `scripts/sources/GSE264489.py` |
| GSE266079 | GEO-associated | GSE266079_genes001 | 315,881 | `scripts/sources/GSE266079.py` |
| GSE266488 | GEO-associated | GSE266488_GSM8178020<br>GSE266488_GSM8178021<br>GSE266488_GSM8178022<br>GSE266488_GSM8178023<br>GSE266488_GSM8178024<br>GSE266488_GSM8178025<br>GSE266488_GSM8178026<br>GSE266488_GSM8178027<br>GSE266488_GSM8178028<br>GSE266488_GSM8178029<br>GSE266488_GSM8178030<br>GSE266488_GSM8248050<br>GSE266488_GSM8248051<br>GSE266488_GSM8248052<br>GSE266488_GSM8248053<br>GSE266488_GSM8248054<br>GSE266488_GSM8248055<br>GSE266488_GSM8248056<br>GSE266488_GSM8248057<br>GSE266488_GSM8248058<br>GSE266488_GSM8248059<br>GSE266488_GSM8248060<br>GSE266488_GSM8248061<br>GSE266488_GSM8248062<br>GSE266488_GSM8248063<br>GSE266488_GSM8248064<br>GSE266488_GSM8248065<br>GSE266488_GSM8248066<br>GSE266488_GSM8248067<br>GSE266488_GSM8248068<br>GSE266488_GSM8248069<br>GSE266488_GSM8248070<br>GSE266488_GSM8248071<br>GSE266488_GSM8248072<br>GSE266488_GSM8248073<br>GSE266488_GSM8248074<br>GSE266488_GSM8248075<br>GSE266488_GSM8248076<br>GSE266488_GSM8248077<br>GSE266488_GSM8248078<br>GSE266488_GSM8248079<br>GSE266488_GSM8248080<br>GSE266488_GSM8248081<br>GSE266488_GSM8248082<br>GSE266488_GSM8248083<br>GSE266488_GSM8248084<br>GSE266488_GSM8248085<br>GSE266488_GSM8248086<br>GSE266488_GSM8248087<br>GSE266488_GSM8248088<br>GSE266488_GSM8248089<br>GSE266488_GSM8248090<br>GSE266488_GSM8248091<br>GSE266488_GSM8248092<br>GSE266488_GSM8248093<br>GSE266488_GSM8248094<br>GSE266488_GSM8248095<br>GSE266488_GSM8248096<br>GSE266488_GSM8248097 | 206,352 | `scripts/sources/GSE266488.py` |
| GSE266519 | GEO-associated | GSE266519 | 88,886 | `scripts/sources/GSE266519.py` |
| GSE267097 | GEO-associated | GSE267097_genes001 | 107,800 | `scripts/sources/GSE267097.py` |
| GSE267644 | GEO-associated | GSE267644 | 100,033 | `scripts/sources/GSE267644.py` |
| GSE267645 | GEO-associated | GSE267645_genes001 | 130,148 | `scripts/sources/GSE267645.py` |
| GSE267718 | GEO-associated | GSE267718_genes001<br>GSE267718_genes002 | 46,761 | `scripts/sources/GSE267718.py` |
| GSE267750 | GEO-associated | GSE267750 | 275,324 | `scripts/sources/GSE267750.py` |
| GSE268210 | GEO-associated | GSE268210_genes001<br>GSE268210_genes002 | 183,738 | `scripts/sources/GSE268210.py` |
| GSE268626 | GEO-associated | GSE268626_genes001 | 26,136 | `scripts/sources/GSE268626.py` |
| GSE268936 | GEO-associated | GSE268936_genes001 | 298,184 | `scripts/sources/GSE268936.py` |
| GSE268995 | GEO-associated | GSE268995_genes001 | 296,318 | `scripts/sources/GSE268995.py` |
| GSE269269 | GEO-associated | GSE269269 | 92,565 | `scripts/sources/GSE269269.py` |
| GSE269338 | GEO-associated | GSE269338 | 25,270 | `scripts/sources/GSE269338.py` |
| GSE269661 | GEO-associated | GSE269661_genes001 | 60,294 | `scripts/sources/GSE269661.py` |
| GSE270165 | GEO-associated | GSE270165 | 53,097 | `scripts/sources/GSE270165.py` |
| GSE270235 | GEO-associated | GSE270235 | 74,032 | `scripts/sources/GSE270235.py` |
| GSE270387 | GEO-associated | GSE270387 | 313,311 | `scripts/sources/GSE270387.py` |
| GSE271284 | GEO-associated | GSE271284_genes001<br>GSE271284_genes002 | 154,940 | `scripts/sources/GSE271284.py` |
| GSE271543 | GEO-associated | GSE271543_genes001 | 35,203 | `scripts/sources/GSE271543.py` |
| GSE271789 | GEO-associated | GSE271789<br>GSE271789_pooled | 447,106 | `scripts/sources/GSE271789.py` |
| GSE272294 | GEO-associated | GSE272294 | 253,519 | `scripts/sources/GSE272294.py` |
| GSE272734 | GEO-associated | GSE272734_genes001<br>GSE272734_genes002<br>GSE272734_genes003<br>GSE272734_genes004<br>GSE272734_genes005<br>GSE272734_genes006<br>GSE272734_genes007<br>GSE272734_genes008<br>GSE272734_genes009<br>GSE272734_genes010<br>GSE272734_genes011<br>GSE272734_genes012<br>GSE272734_genes013 | 160,684 | `scripts/sources/GSE272734.py` |
| GSE272735 | GEO-associated | GSE272735_genes001<br>GSE272735_genes002<br>GSE272735_genes003<br>GSE272735_genes004 | 30,798 | `scripts/sources/GSE272735.py` |
| GSE272993 | GEO-associated | GSE272993_additional_selected_cd8<br>GSE272993_pbmc_selected_cd8 | 117,970 | `scripts/sources/GSE272993.py` |
| GSE273718 | GEO-associated | GSE273718_genes001 | 17,838 | `scripts/sources/GSE273718.py` |
| GSE273720 | GEO-associated | GSE273720_genes001 | 50,229 | `scripts/sources/GSE273720.py` |
| GSE273781 | GEO-associated | GSE273781_genes001 | 115,903 | `scripts/sources/GSE273781.py` |
| GSE274284 | GEO-associated | GSE274284 | 114,707 | `scripts/sources/GSE274284.py` |
| GSE274381 | GEO-associated | GSE274381_genes001<br>GSE274381_genes002<br>GSE274381_genes003<br>GSE274381_genes004 | 54,730 | `scripts/sources/GSE274381.py` |
| GSE274588 | GEO-associated | GSE274588 | 219,270 | `scripts/sources/GSE274588.py` |
| GSE275071 | GEO-associated | GSE275071 | 17,350 | `scripts/sources/GSE275071.py` |
| GSE275259 | GEO-associated | GSE275259_genes001 | 62,419 | `scripts/sources/GSE275259.py` |
| GSE275947 | GEO-associated | GSE275947_genes001 | 135,234 | `scripts/sources/GSE275947.py` |
| GSE276742 | GEO-associated | GSE276742 | 88,900 | `scripts/sources/GSE276742.py` |
| GSE277117 | GEO-associated | GSE277117 | 311,993 | `scripts/sources/GSE277117.py` |
| GSE277496 | GEO-associated | GSE277496_RNA | 10,489 | `scripts/sources/GSE277496.py` |
| GSE277590 | GEO-associated | GSE277590 | 162,004 | `scripts/sources/GSE277590.py` |
| GSE278119 | GEO-associated | GSE278119 | 216,513 | `scripts/sources/GSE278119.py` |
| GSE278838 | GEO-associated | GSE278838_genes001 | 279,317 | `scripts/sources/GSE278838.py` |
| GSE279452 | GEO-associated | GSE279452 | 942,199 | `scripts/sources/GSE279452.py` |
| GSE279838 | GEO-associated | GSE279838_genes001 | 28,080 | `scripts/sources/GSE279838.py` |
| GSE279945 | GEO-associated | GSE279945_author_sc_counts_reannotated_with_counts | 298,083 | `scripts/sources/GSE279945.py` |
| GSE280258 | GEO-associated | GSE280258_genes001 | 5,164 | `scripts/sources/GSE280258.py` |
| GSE280401 | GEO-associated | GSE280401_genes001 | 7,595 | `scripts/sources/GSE280401.py` |
| GSE280720 | GEO-associated | GSE280720_genes001 | 182,717 | `scripts/sources/GSE280720.py` |
| GSE280982 | GEO-associated | GSE280982_genes001 | 34,246 | `scripts/sources/GSE280982.py` |
| GSE281237 | GEO-associated | GSE281237_genes001 | 236,143 | `scripts/sources/GSE281237.py` |
| GSE282427 | GEO-associated | GSE282427_genes001 | 163,148 | `scripts/sources/GSE282427.py` |
| GSE282793 | GEO-associated | GSE282793_genes001 | 44,413 | `scripts/sources/GSE282793.py` |
| GSE283397 | GEO-associated | GSE283397 | 35,699 | `scripts/sources/GSE283397.py` |
| GSE283746 | GEO-associated | GSE206289_pCoV1<br>GSE283746 | 440,148 | `scripts/sources/GSE283746.py` |
| GSE284075 | GEO-associated | GSE284075_genes001 | 72,970 | `scripts/sources/GSE284075.py` |
| GSE284677 | GEO-associated | GSE284677_peripheral_genes001 | 134,515 | `scripts/sources/GSE284677.py` |
| GSE285173 | GEO-associated | GSE285173_GSM8697592<br>GSE285173_GSM8697593<br>GSE285173_GSM8697594<br>GSE285173_GSM8697595<br>GSE285173_GSM8697596<br>GSE285173_GSM8697597 | 75,332 | `scripts/sources/GSE285173.py` |
| GSE285335 | GEO-associated | GSE285335_genes001 | 281,513 | `scripts/sources/GSE285335.py` |
| GSE285646 | GEO-associated | GSE285646_genes001 | 156,089 | `scripts/sources/GSE285646.py` |
| GSE285773 | GEO-associated | GSE285773_genes001<br>GSE285773_genes002 | 262,929 | `scripts/sources/GSE285773.py` |
| GSE285774 | GEO-associated | GSE285774_genes001 | 46,233 | `scripts/sources/GSE285774.py` |
| GSE286325 | GEO-associated | GSE286325 | 55,469 | `scripts/sources/GSE286325.py` |
| GSE287541 | GEO-associated | GSE287541 | 173,082 | `scripts/sources/GSE287541.py` |
| GSE288147 | GEO-associated | GSE288147_genes001 | 161,318 | `scripts/sources/GSE288147.py` |
| GSE288581 | GEO-associated | GSE288581_genes001 | 17,824 | `scripts/sources/GSE288581.py` |
| GSE288613 | GEO-associated | GSE288613_genes001 | 126,887 | `scripts/sources/GSE288613.py` |
| GSE290026 | GEO-associated | GSE290026 | 110,055 | `scripts/sources/GSE290026.py` |
| GSE290061 | GEO-associated | GSE290061_genes001 | 31,831 | `scripts/sources/GSE290061.py` |
| GSE290557 | GEO-associated | GSE290557_genes001 | 99,933 | `scripts/sources/GSE290557.py` |
| GSE290679 | GEO-associated | GSE290679 | 644,145 | `scripts/sources/GSE290679.py` |
| GSE290722 | GEO-associated | GSE290722 | 405,775 | `scripts/sources/GSE290722.py` |
| GSE290836 | GEO-associated | GSE290836 | 351,299 | `scripts/sources/GSE290836.py` |
| GSE291111 | GEO-associated | GSE291111 | 4,326 | `scripts/sources/GSE291111.py` |
| GSE292141 | GEO-associated | GSE292141_genes001 | 130,165 | `scripts/sources/GSE292141.py` |
| GSE293860 | GEO-associated | GSE293860_genes001 | 123,483 | `scripts/sources/GSE293860.py` |
| GSE293960 | GEO-associated | GSE293960_genes001 | 126,777 | `scripts/sources/GSE293960.py` |
| GSE294268 | GEO-associated | GSE294268_genes001 | 142,054 | `scripts/sources/GSE294268.py` |
| GSE295129 | GEO-associated | GSE295129_genes001 | 31,152 | `scripts/sources/GSE295129.py` |
| GSE295173 | GEO-associated | GSE295173 | 13,879 | `scripts/sources/GSE295173.py` |
| GSE295174 | GEO-associated | GSE295174 | 2,299 | `scripts/sources/GSE295174.py` |
| GSE295359 | GEO-associated | GSE295359 | 478,079 | `scripts/sources/GSE295359.py` |
| GSE295491 | GEO-associated | GSE295491_genes001 | 309,423 | `scripts/sources/GSE295491.py` |
| GSE295527 | GEO-associated | GSE295527_archive_genes001<br>GSE295527_archive_genes002<br>GSE295527_flex_gsm8951957<br>GSE295527_flex_gsm8951958<br>GSE295527_flex_gsm8951959<br>GSE295527_flex_gsm8951960<br>GSE295527_genes001 | 625,888 | `scripts/sources/GSE295527.py` |
| GSE295601 | GEO-associated | GSE295601_genes001 | 23,561 | `scripts/sources/GSE295601.py` |
| GSE296678 | GEO-associated | GSE296678_genes001 | 46,128 | `scripts/sources/GSE296678.py` |
| GSE296867 | GEO-associated | GSE296867_genes001 | 45,607 | `scripts/sources/GSE296867.py` |
| GSE297529 | GEO-associated | GSE297529_genes001 | 43,069 | `scripts/sources/GSE297529.py` |
| GSE297547 | GEO-associated | GSE297547_C001<br>GSE297547_C002<br>GSE297547_C003<br>GSE297547_C004<br>GSE297547_C005<br>GSE297547_C006<br>GSE297547_C007<br>GSE297547_C008<br>GSE297547_C009<br>GSE297547_C010<br>GSE297547_MG001<br>GSE297547_MG002<br>GSE297547_MG003<br>GSE297547_MG004<br>GSE297547_MG005<br>GSE297547_MG006<br>GSE297547_MG007<br>GSE297547_MG008<br>GSE297547_MG009 | 242,232 | `scripts/sources/GSE297547.py` |
| GSE298838 | GEO-associated | GSE298838 | 159,495 | `scripts/sources/GSE298838.py` |
| GSE300837 | GEO-associated | GSE300837 | 97,157 | `scripts/sources/GSE300837.py` |
| GSE301651 | GEO-associated | GSE301651_genes001 | 8,351 | `scripts/sources/GSE301651.py` |
| GSE302810 | GEO-associated | GSE302810_genes001 | 24,551 | `scripts/sources/GSE302810.py` |
| GSE303333 | GEO-associated | GSE303333_genes001 | 126,899 | `scripts/sources/GSE303333.py` |
| GSE303537 | GEO-associated | GSE303537_genes001 | 575,112 | `scripts/sources/GSE303537.py` |
| GSE304309 | GEO-associated | GSE304309 | 75,147 | `scripts/sources/GSE304309.py` |
| GSE304657 | GEO-associated | GSE304657_genes001 | 540,549 | `scripts/sources/GSE304657.py` |
| GSE304704 | GEO-associated | GSE304704_genes001 | 92,508 | `scripts/sources/GSE304704.py` |
| GSE304878 | GEO-associated | GSE304878_genes001 | 98,650 | `scripts/sources/GSE304878.py` |
| GSE304891 | GEO-associated | GSE304891_genes001 | 232,592 | `scripts/sources/GSE304891.py` |
| GSE305477 | GEO-associated | GSE305477 | 3,385 | `scripts/sources/GSE305477.py` |
| GSE305637 | GEO-associated | GSE305637_genes001 | 57,178 | `scripts/sources/GSE305637.py` |
| GSE306003 | GEO-associated | GSE306003_genes001 | 47,389 | `scripts/sources/GSE306003.py` |
| GSE306378 | GEO-associated | GSE306378_genes001 | 78,127 | `scripts/sources/GSE306378.py` |
| GSE306441 | GEO-associated | GSE306441 | 124,173 | `scripts/sources/GSE306441.py` |
| GSE308075 | GEO-associated | GSE308075_genes001 | 113,509 | `scripts/sources/GSE308075.py` |
| GSE308745 | GEO-associated | GSE308745_genes001 | 224,011 | `scripts/sources/GSE308745.py` |
| GSE309140 | GEO-associated | GSE309140 | 27,868 | `scripts/sources/GSE309140.py` |
| GSE309247 | GEO-associated | GSE309247 | 139,122 | `scripts/sources/GSE309247.py` |
| GSE309368 | GEO-associated | GSE309368_genes001 | 79,406 | `scripts/sources/GSE309368.py` |
| GSE309595 | GEO-associated | GSE309595_genes001 | 4,019,103 | `scripts/sources/GSE309595.py` |
| GSE309906 | GEO-associated | GSE309906_genes001 | 116,624 | `scripts/sources/GSE309906.py` |
| GSE309970 | GEO-associated | GSE309970_final | 170,907 | `scripts/sources/GSE309970.py` |
| GSE310662 | GEO-associated | GSE310662 | 105,181 | `scripts/sources/GSE310662.py` |
| GSE310935 | GEO-associated | GSE310935_genes001 | 31,764 | `scripts/sources/GSE310935.py` |
| GSE311112 | GEO-associated | GSE311112_genes001 | 140,082 | `scripts/sources/GSE311112.py` |
| GSE312336 | GEO-associated | GSE312336_genes001 | 24,304 | `scripts/sources/GSE312336.py` |
| GSE313043 | GEO-associated | GSE313043_genes001<br>GSE313043_pooled | 108,252 | `scripts/sources/GSE313043.py` |
| GSE313642 | GEO-associated | GSE313642_genes001 | 239,347 | `scripts/sources/GSE313642.py` |
| GSE313894 | GEO-associated | GSE313894_genes001 | 86,898 | `scripts/sources/GSE313894.py` |
| GSE314010 | GEO-associated | GSE314010_genes001 | 52,820 | `scripts/sources/GSE314010.py` |
| GSE314416 | GEO-associated | GSE314416_genes001 | 566,024 | `scripts/sources/GSE314416.py` |
| GSE314863 | GEO-associated | GSE314863_5prime | 17,574 | `scripts/sources/GSE314863.py` |
| GSE316301 | GEO-associated | GSE316301_GSM9449215<br>GSE316301_GSM9449216<br>GSE316301_GSM9449217<br>GSE316301_GSM9449218<br>GSE316301_GSM9449219<br>GSE316301_GSM9449220<br>GSE316301_GSM9449221<br>GSE316301_GSM9449222<br>GSE316301_GSM9449223<br>GSE316301_GSM9449224<br>GSE316301_GSM9449225<br>GSE316301_GSM9449226<br>GSE316301_GSM9449227<br>GSE316301_GSM9449228<br>GSE316301_GSM9449229<br>GSE316301_GSM9449230<br>GSE316301_GSM9449231<br>GSE316301_GSM9449232<br>GSE316301_GSM9449233<br>GSE316301_GSM9449234<br>GSE316301_GSM9449235<br>GSE316301_GSM9449236<br>GSE316301_GSM9449237<br>GSE316301_GSM9449238 | 309,713 | `scripts/sources/GSE316301.py` |
| GSE316337 | GEO-associated | GSE316337 | 376,996 | `scripts/sources/GSE316337.py` |
| GSE316688 | GEO-associated | GSE316688_discovery<br>GSE316688_nonresponder | 288,707 | `scripts/sources/GSE316688.py` |
| GSE317309 | GEO-associated | GSE317309_d0_genes001<br>GSE317309_d21_genes001<br>GSE317309_d42_genes001<br>GSE317309_d63_genes001 | 258,254 | `scripts/sources/GSE317309.py` |
| GSE317605 | GEO-associated | GSE317605_genes001 | 52,388 | `scripts/sources/GSE317605.py` |
| GSE317902 | GEO-associated | GSE317902 | 72,867 | `scripts/sources/GSE317902.py` |
| GSE318371 | GEO-associated | GSE318371 | 376,051 | `scripts/sources/GSE318371.py` |
| GSE318946 | GEO-associated | GSE318946_genes001 | 173,895 | `scripts/sources/GSE318946.py` |
| GSE319710 | GEO-associated | GSE319710 | 73,865 | `scripts/sources/GSE319710.py` |
| GSE320359 | GEO-associated | GSE320359_genes001<br>GSE320359_genes002 | 18,817 | `scripts/sources/GSE320359.py` |
| GSE324815 | GEO-associated | GSE324815_genes001 | 50,817 | `scripts/sources/GSE324815.py` |
| GSE324824 | GEO-associated | GSE324824 | 79,393 | `scripts/sources/GSE324824.py` |
| GSE326288 | GEO-associated | GSE326288_genes001 | 310,571 | `scripts/sources/GSE326288.py` |
| GSE326661 | GEO-associated | GSE326661_genes001 | 92,102 | `scripts/sources/GSE326661.py` |
| GSE328561 | GEO-associated | GSE328561_genes001<br>GSE328561_genes002 | 78,178 | `scripts/sources/GSE328561.py` |
| GSE328788 | GEO-associated | GSE328788_genes001<br>GSE328788_genes002 | 108,862 | `scripts/sources/GSE328788.py` |
| GSE329190 | GEO-associated | GSE329190_genes001 | 72,067 | `scripts/sources/GSE329190.py` |
| GSE330576 | GEO-associated | GSE330576_genes001 | 99,456 | `scripts/sources/GSE330576.py` |
| GSE332711 | GEO-associated | GSE332711_genes001 | 175,746 | `scripts/sources/GSE332711.py` |
| GSE332837 | GEO-associated | GSE332837_genes001 | 273,269 | `scripts/sources/GSE332837.py` |
| GSE334503 | GEO-associated | GSE334503_author_native | 322,174 | `scripts/sources/GSE334503.py` |
| GSE334857 | GEO-associated | GSE334857_conditional_views | 12,835 | `scripts/sources/GSE334857.py` |
| GSE335494 | GEO-associated | GSE335494_genes001<br>GSE335494_genes002 | 204,982 | `scripts/sources/GSE335494.py` |
| GSE335977 | GEO-associated | GSE335977_genes001 | 68,692 | `scripts/sources/GSE335977.py` |
| GSE335997 | GEO-associated | GSE335997_genes001 | 34,841 | `scripts/sources/GSE335997.py` |
| GSE336230 | GEO-associated | GSE336230 | 32,752 | `scripts/sources/GSE336230.py` |
| GSE337705 | GEO-associated | GSE337705_genes001<br>GSE337705_genes002<br>GSE337705_genes003 | 266,603 | `scripts/sources/GSE337705.py` |
| GSE338165 | GEO-associated | GSE338165_genes001 | 264,972 | `scripts/sources/GSE338165.py` |
| GSE339254 | GEO-associated | GSE339254_genes001 | 174,703 | `scripts/sources/GSE339254.py` |
| GSE339599 | GEO-associated | GSE339599_genes001 | 16,771 | `scripts/sources/GSE339599.py` |
| GSE341102 | GEO-associated | GSE341102_genes001 | 68,976 | `scripts/sources/GSE341102.py` |
| GSE341191 | GEO-associated | GSE341191_genes001 | 198,150 | `scripts/sources/GSE341191.py` |
| GSE342939 | GEO-associated | GSE342939_genes001 | 195,168 | `scripts/sources/GSE342939.py` |
| GSE345406 | GEO-associated | GSE345406_genes001 | 43,189 | `scripts/sources/GSE345406.py` |
| GSE345981 | GEO-associated | GSE345981 | 66,650 | `scripts/sources/GSE345981.py` |
| GSE347367 | GEO-associated | GSE347367_BNT_V07_ID0<br>GSE347367_BNT_V07_ID1<br>GSE347367_BNT_V09_ID0<br>GSE347367_BNT_V09_ID1<br>GSE347367_BNT_V30_ID0<br>GSE347367_BNT_V30_ID1<br>GSE347367_BNT_V30_ID2<br>GSE347367_BNT_V43_ID0<br>GSE347367_BNT_V43_ID1<br>GSE347367_BNT_V43_ID3<br>GSE347367_BNT_V46_ID0<br>GSE347367_BNT_V46_ID1<br>GSE347367_BNT_V46_ID3<br>GSE347367_BNT_V69_ID0<br>GSE347367_BNT_V69_ID1<br>GSE347367_BNT_V69_ID2<br>GSE347367_BNT_V69_ID3<br>GSE347367_BNT_V74_ID0<br>GSE347367_BNT_V74_ID1<br>GSE347367_BNT_V74_ID3<br>GSE347367_Moderna_V01_ID1<br>GSE347367_Moderna_V02_ID0<br>GSE347367_Moderna_V10_ID1<br>GSE347367_Moderna_V10_ID2<br>GSE347367_Moderna_V10_ID3<br>GSE347367_Moderna_V14_ID1<br>GSE347367_Moderna_V16_ID0<br>GSE347367_Moderna_V16_ID1<br>GSE347367_Moderna_V20_ID0<br>GSE347367_Moderna_V20_ID1<br>GSE347367_Moderna_V23_ID0<br>GSE347367_Moderna_V23_ID1<br>GSE347367_Moderna_V23_ID3<br>GSE347367_Moderna_V24_ID0<br>GSE347367_Moderna_V24_ID1<br>GSE347367_Moderna_V26_ID0<br>GSE347367_Moderna_V26_ID1<br>GSE347367_Moderna_V26_ID3<br>GSE347367_Moderna_V27_ID1<br>GSE347367_Moderna_V36_ID0<br>GSE347367_Moderna_V36_ID1<br>GSE347367_Moderna_V38_ID1<br>GSE347367_Moderna_V41_ID1<br>GSE347367_Moderna_V44_ID0<br>GSE347367_Moderna_V44_ID1<br>GSE347367_Moderna_V44_ID3<br>GSE347367_Moderna_V45_ID0<br>GSE347367_Moderna_V45_ID1<br>GSE347367_Moderna_V45_ID3<br>GSE347367_Moderna_V47_ID0<br>GSE347367_Moderna_V47_ID1<br>GSE347367_Moderna_V47_ID2<br>GSE347367_Moderna_V50_ID0<br>GSE347367_Moderna_V50_ID1<br>GSE347367_Moderna_V50_ID2<br>GSE347367_Moderna_V51_ID0<br>GSE347367_Moderna_V52_ID0<br>GSE347367_Moderna_V52_ID1<br>GSE347367_Moderna_V53_ID0<br>GSE347367_Moderna_V53_ID1<br>GSE347367_Moderna_V53_ID3<br>GSE347367_Moderna_V54_ID0<br>GSE347367_Moderna_V54_ID1<br>GSE347367_Moderna_V54_ID2<br>GSE347367_Moderna_V54_ID3<br>GSE347367_Moderna_V55_ID1<br>GSE347367_Moderna_V56_ID0<br>GSE347367_Moderna_V56_ID1<br>GSE347367_Moderna_V57_ID0<br>GSE347367_Moderna_V57_ID1<br>GSE347367_Moderna_V60_ID1<br>GSE347367_Moderna_V60_ID2<br>GSE347367_Moderna_V60_ID3<br>GSE347367_Moderna_V61_ID1<br>GSE347367_Moderna_V61_ID3<br>GSE347367_Moderna_V65_ID1<br>GSE347367_Moderna_V77_ID1 | 190,101 | `scripts/sources/GSE347367.py` |
| GSE347634 | GEO-associated | GSE347634 | 281,799 | `scripts/sources/GSE347634.py` |
| GSE349910 | GEO-associated | GSE349910_genes001 | 243,853 | `scripts/sources/GSE349910.py` |
| GSE350157 | GEO-associated | GSE350157 | 40,822 | `scripts/sources/GSE350157.py` |
| KCL_immune_aging2026 | Author-hosted atlas resources | KCL_immune_aging2026 | 758,026 | `scripts/sources/KCL_immune_aging2026.py` |
| PRJEB86440 | ENA | PRJEB86440 | 309,564 | `scripts/sources/PRJEB86440.py` |
| SPAC_immune_aging2026 | Author-hosted atlas resources | SPAC_immune_aging2026 | 293,707 | `scripts/sources/SPAC_immune_aging2026.py` |
| c672834e-c3e3-49cb-81a5-4c844be4a975 | CELLxGENE | a199ca73-035d-44e2-9893-4c493151db21 | 105,785 | `scripts/sources/c672834e-c3e3-49cb-81a5-4c844be4a975.py` |
| cc431242-35ea-41e1-a100-41e0dec2665b | CELLxGENE | 1b350d0a-4535-4879-beb6-1142f3f94947 | 111,319 | `scripts/sources/cc431242-35ea-41e1-a100-41e0dec2665b.py` |
| ced320a1-29f3-47c1-a735-513c7084d508 | CELLxGENE | c838aec3-03ef-4398-b882-0e3912abfff0 | 1,265,167 | `scripts/sources/ced320a1-29f3-47c1-a735-513c7084d508.py` |
| d1e0e64d-6d2a-4a3e-b7f4-43ed909a9d9c | CELLxGENE | ca7d95ac-53aa-4a94-8c76-300c608ce6a0 | 435,176 | `scripts/sources/d1e0e64d-6d2a-4a3e-b7f4-43ed909a9d9c.py` |
| dde06e0f-ab3b-46be-96a2-a8082383c4a1 | CELLxGENE | 3faad104-2ab8-4434-816d-474d8d2641db | 1,248,980 | `scripts/sources/dde06e0f-ab3b-46be-96a2-a8082383c4a1.py` |
| ddfad306-714d-4cc0-9985-d9072820c530 | CELLxGENE | c7775e88-49bf-4ba2-a03b-93f00447c958 | 640,901 | `scripts/sources/ddfad306-714d-4cc0-9985-d9072820c530.py` |
| e1a9ca56-f2ee-435d-980a-4f49ab7a952b | CELLxGENE | d18736c3-6292-4379-919a-d6d973204c87 | 45,950 | `scripts/sources/e1a9ca56-f2ee-435d-980a-4f49ab7a952b.py` |

## Analysis conventions

- Missing cell types are never replaced by zero-expression profiles.
- Duplicate gene symbols are averaged after conversion to symbols.
- Genes are selected by summed within-map variance, before mean imputation.
- Each reference profile is the unweighted mean of its observed stars.
- Alignment uses one equal-star translation per map and is not target-aware.
- UMAP is visualization only; Figure 2 is evaluated in 2,000-gene space.
