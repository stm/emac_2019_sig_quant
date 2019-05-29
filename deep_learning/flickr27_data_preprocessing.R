# data preprocessing
# 
# - download the file set from here: http://image.ntua.gr/iva/datasets/flickr_logos/
# - extract image to a folder (e.g., flickr27)
# -copy files into the right directories

rm(list=ls())

# make results reproducible
set.seed(2787)

library(dplyr)

# read files
options(stringsAsFactors = F)
df_train <- read.table("flickr27/flickr_logos_27_dataset_training_set_annotation.txt", sep=" ")
df_test <- read.table("flickr27/flickr_logos_27_dataset_query_set_annotation.txt", sep="\t")
names(df_test) <- c("filename", "label")

train_files <- list.files("flickr27/flickr_logos_27_dataset_images_cropped/")
df_train$filename <- paste0(1:nrow(df_train),"_",df_train$V1)

# check that  training files are not part of test set
stopifnot(all(!(df_test$V1 %in% train_files)))
# check whether cropped images are all in original data
stopifnot(all(df_train$filename %in% train_files))

# only keep relevant information in training set
df_train <- df_train[c(1,2,9)]
names(df_train) <- c("filename_orig", "label", "filename")

# how many images per class?
table(df_train$label)

# only keep Coca Cola, Pepsi, and Heineken in data sets
df_test <- df_test %>% filter(label == "Cocacola" | label == "Pepsi" | label == "Heineken")
df_train <- df_train %>% filter(label == "Cocacola" | label == "Pepsi" | label == "Heineken")

# do 70:30 split on training data to create validaion data
train <- df_train %>% group_by(label) %>% sample_n(120)
val <- df_train[!(df_train$filename %in% train$filename),]
val <- val %>% group_by(label) %>% sample_n(40)

# check whether validation and training data do not overlap
stopifnot(all(!(train$filename %in% val$filename)))


# copy files
for(i in unique(train$label)){
  dir.create(paste0("flickr27/training/", i), recursive = TRUE)
}

for( i in unique(val$label)){
  dir.create(paste0("flickr27/validation/", i), recursive = TRUE)
}

for( i in unique(df_test$label)){
  dir.create(paste0("flickr27/test/", i), recursive = TRUE)
}


# test set
for(i in df_test$filename){
  print(i)
  label <- df_test[df_test$filename == i,]$label
  file.copy(paste0("flickr27/flickr_logos_27_dataset_images/", i), paste0("flickr27/test/", label, "/", i))
}

# validation set
for(i in val$filename){
  print(i)
  label <- val[val$filename == i,]$label
  file.copy(paste0("flickr27/flickr_logos_27_dataset_images_cropped/", i), paste0("flickr27/validation/", label, "/", i))
}

# validation set
for(i in train$filename){
  print(i)
  label <- train[train$filename == i,]$label
  file.copy(paste0("flickr27/flickr_logos_27_dataset_images_cropped/", i), paste0("flickr27/training/", label, "/", i))
}
