# Coke vs. Pepsi: A deep learning brand logo classifier using transfer learning
#
# use transfer learning to train a convolutional neural network for binary image
# classification (pepsi vs coke brand logo detection)
#
# Stefan Mayer, University of Tuebingen, May 2019

# clean workspace
rm(list=ls())

# load keras library
library(keras)

# make results reproducible
myseed <- 42
set.seed(myseed)
use_session_with_seed(myseed, disable_gpu = FALSE, disable_parallel_cpu = FALSE)

# clear keras session (if neccessary)
k_clear_session()

# set data directories (i.e., location of images)
#
# beware of the folder structure!
# the labels (i.e., brand names) will be inferred from the folder names
train_directory <- "flickr27/training/"
val_directory <- "flickr27/validation/"
test_directory <- "flickr27/test/"

# how many images in training, validatiton and test set
train_samples = length(list.files(path = train_directory, recursive = TRUE))
validation_samples = length(list.files(path = val_directory, recursive = TRUE))
test_samples = length(list.files(path = test_directory, recursive = TRUE))

# set parameters: image width and height to use for the model
img_width <- 100
img_height <- 100

# batch size: how many samples (i.e., images) should be used in one training iteration
batch_size <- 10

# we will sample images from the folder (instead of loading all images at once)
#
# we can define 'image_data_generator' function that performes operations when
# loading the images (here, scale the values between 0 and 1)
datagen_train <- image_data_generator(rescale = 1/255)
datagen_val <- image_data_generator(rescale = 1/255)

# in this example, we'll only use the data from Pepsi and Coke
class_list <- c("Pepsi", "Cocacola")

# the 'flow_images_from_directory' function samples images from a folder with
# the given parameters (e.g., data generator function, image size, batch size, ...)
# 
# note that we use "binary" as class mode b/c we have a binary classification task
train_generator <- flow_images_from_directory(train_directory, generator = datagen_train,
                                              target_size = c(img_width, img_height),
                                              class_mode = "binary", batch_size = batch_size,
                                              classes = class_list,
                                              seed = myseed)

validation_generator <- flow_images_from_directory(val_directory, generator = datagen_val,
                                                   target_size = c(img_width, img_height),
                                                   class_mode = "binary", batch_size = batch_size,
                                                   classes = class_list,
                                                   seed = myseed)

# sanity check: were all images detected?
# 
# check label coding
train_generator$class_indices
# 
# note that coke is coded as 1 and pepsi as 0
table(train_generator$classes)



# define our neural network
# 
# as we use transfer learning, we have to decide on the base architecture on
# top of which we'll build our new classifier
#
# here, we use the popular VGG16 network pre-trained on imagenet
# luckily, this comes pre-packaged with keras.
#
# VGG16 as convolutional base with the imagenet weights
conv_base <- application_vgg16(
  weights = "imagenet",
  include_top = FALSE,
  input_shape = c(img_width, img_height, 3)
)

# inspect model
summary(conv_base)

# now we add our custom layers (i.e., a new classifier on top of the
# convolutional base)
model <- keras_model_sequential() %>%
  conv_base %>%
  layer_flatten() %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")

# freeze the convolutional base (otherwise, we would re-train the whole model)
cat("Trainable weights before freezing:",length(model$trainable_weights), "\n")

freeze_weights(conv_base)

cat("Trainable weights after freezing:",length(model$trainable_weights), "\n")


# visualize the model
# 
# note that you need to have deepviz intalled
# you can get it e.g. from github
# devtools::install_github("andrie/deepviz")
deepviz::plot_model(model)


# compile the model
# 
# note that we use "binary_crossentropy" as loss function because of our binary
# classification task
model %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_rmsprop(lr = 5e-5, decay = 1e-6),
  metrics = "accuracy"
)

# train the model (this may take several minutes ...)
# 
# note that the number of epochs defines how many iterations should be done.
hist <- model %>% fit_generator(
  train_generator,
  steps_per_epoch = as.integer(train_samples/batch_size), 
  epochs = 10,
  validation_data = validation_generator,
  validation_steps = as.integer(validation_samples/batch_size)
)



# now that our models is trained, we can evaluate the model's performance
# 
# first, let's have a look at two generic images downloaded from google. we
# define a image prediction function that does some image preprocessing (so that
# the data fit the data format with which the network was trained) and predicts
# the probability that the image contains a coke logo (coke was coded as 1 and
# pepsi as zero, see above)
pred_img <- function(path) {
  img <- image_load(path, target_size = c(img_width, img_height))
  x <- image_to_array(img)
  x <- x/255
  x <- array_reshape(x, c(1, dim(x)))
  return(paste0(round(model %>% predict(x)*100,3),"%"))
}

# also create a helper function to plot an image
img_plot <- function(path){
  img <- image_load(path, target_size = c(img_width, img_height))
  x <- image_to_array(img)
  x <- x/255
  grid::grid.raster(x)
}

# reminder: Coca Cola: 1, Pepsi: 0
# pred_img gives prob of the image containing a coca Cola logo
img_plot("cocacola.jpg")
pred_img("cocacola.jpg")

img_plot("pepsi.jpg")
pred_img("pepsi.jpg")

# performance on official test set
test_generator <- flow_images_from_directory(
  test_directory, generator = datagen_val,
  target_size = c(img_width, img_height),
  class_mode = "binary", batch_size = batch_size,
  classes = class_list,
  seed = myseed)

test_performance <- model %>% evaluate_generator(test_generator, steps = 5)
print(paste0("Test accuracy: ", round(test_performance$acc*100,4), "%"))

