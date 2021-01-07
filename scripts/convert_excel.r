library(openxlsx)
# The ```assay``` function allows us to covert deseq2 objects into data tables
data = assay(se)
# The ```write.xlsx``` function allows to convert datatables into to excel files
write.xlsx(data, file='se.xlsx', sheetName="Sheet1",
           col.names=TRUE, row.names=TRUE, append=FALSE)