# The following code is meant to set up the connection to the AWS S3
# https://www.gormanalysis.com/blog/connecting-to-aws-s3-with-r/
# Note that the guide above creates a new user specifically for RStudio. For our purposes, we won't do this unless we are hosting this for someone outside the lab. 

# https://github.com/cloudyr/aws.s3
install.packages("aws.s3")
library(aws.s3)

# PLEASE NOTE that the RStudio Server WILL NOT save your R Environment if you STOP this instance. Please save all in-memory R objects AT LEAST to local storage and PREFERABLY to S3 before shutting down.
# Local storage will be preserved if you stop the instance.

# Below are the credentials for a user account that ONLY has access to the lab-tp-rstudio-scratch bucket

#   "AWS_SECRET_ACCESS_KEY" =
#   "AWS_DEFAULT_REGION" = "us-east-1")

# List of all available S3 buckets. Note that while you are able to list all buckets, you will only have access to lab-tp-rstudio-scratch
bucketlist()

# Get a listing of all objects in a bucket
get_bucket(bucket = "lab-tp-rstudio-scratch/palomerolab/R")

# To upload files - save files locally, and then use put_object() to upload to the bucket of interest
# Please use your own folder within the lab-tp-rstudio-scratch bucket as shown below.
put_object(
  file = "/home/palomerolab/temp.txt", 
  object = "temp.txt", 
  bucket = "lab-tp-rstudio-scratch/palomerolab"
)

# There are two methods of getting objects from S3
# (1) Download objects locally and then read them into R using whatever functions you may normally use
save_object("temp.txt", file = "temp.txt", bucket = "lab-tp-rstudio-scratch/palomerolab")

# (2) Read objects directly into R memory as you normally would using s3read_using()
temp <- s3read_using(FUN = read.table, # The function you would normally use to read in object from local
             object = "temp.txt",
             bucket = "lab-tp-rstudio-scratch/palomerolab")

# Using AWS CLI is recommended for syncing multiple files at once. This can only be done in Terminal, which is available in R Studio as a tab next to "Console"
# $ aws s3 sync s3://lab-tp-rstudio-scratch/palomerolab/ /home/palomerolab

# The aws.s3 package also has specific functions for in-memory R objects
# s3save is analogous to save()
s3save(
  RObject_Name,
  object = "R_Object_Name.Rdata",
  bucket = "lab-tp-rstudio-scratch/palomerolab"
)

# s3saveRDS is analogous to saveRDS()
s3save(
  RObject_Name,
  object = "R_Object_Name.rds",
  bucket = "lab-tp-rstudio-scratch/palomerolab"
)

# s3load is analogous to load()
s3load(
  RObject_Name,
  object = "R_Object_Name.Rdata",
  bucket = "lab-tp-rstudio-scratch/palomerolab"
)

# s3readRDS is analogous to readRDS()
s3readRDS(
  RObject_Name,
  object = "R_Object_Name.rds",
  bucket = "lab-tp-rstudio-scratch/palomerolab"
)
