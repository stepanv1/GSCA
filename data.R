#creation of data objects to test functions 
W0<-read.table("./temp/W_recoded_data.txt")
B0<-read.table("./temp/B.txt")
Z0<-read.table("./temp/Data_recoded_Z.txt")

save(list = ls(all=TRUE), file= "test.Rdata")
