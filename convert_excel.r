library(openxlsx)
# The ```assay``` function allows us to covert deseq2 objects into data tables
data = assay(se)
# The ```write.xlsx``` function allows to convert datatables into to excel files
write.xlsx(data,rownames = TRUE, file='example.xlsx')