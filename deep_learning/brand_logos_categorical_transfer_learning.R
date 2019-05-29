# Coke vs. Pepsi vs. Heineken
# A deep learning brand logo classifier using transfer learning
#
# use transfer learning to train a convolutional neural network for multiclass
# image classification (pepsi vs. coke vs. heineken brand logo detection)
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

# the 'flow_images_from_directory' function samples images from a folder with
# the given parameters (e.g., data generator function, image size, batch size, ...)
# 
# note that we use "categorical" as class mode b/c we have a multiclass classification task
train_generator <- flow_images_from_directory(train_directory, generator = datagen_train,
                                              target_size = c(img_width, img_height),
                                              class_mode = "categorical", batch_size = batch_size,
                                              seed = myseed)

validation_generator <- flow_images_from_directory(val_directory, generator = datagen_val,
                                                   target_size = c(img_width, img_height),
                                                   class_mode = "categorical", batch_size = batch_size,
                                                   seed = myseed)

# sanity check: were all images detected?
# 
# check label coding
train_generator$class_indices
#
# store number of classes in variable
num_classes <- length(train_generator$class_indices)
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
  layer_dense(units = 256, activation = "relu") %>%
  layer_dense(units = num_classes, activation = "softmax")

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
# note that we use "categorical_crossentropy" as loss function
# because of our multiclass classification task
model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = optimizer_rmsprop(lr = 1e-5, decay = 1e-6),
  metrics = "accuracy"
)

# train the model (this may take some time ...)
# 
# note that the number of epochs defines how many iterations should be done.
hist <- model %>% fit_generator(
  train_generator,
  steps_per_epoch = as.integer(train_samples/batch_size), 
  epochs = 20,
  validation_data = validation_generator,
  validation_steps = as.integer(validation_samples/batch_size)
)



# now that our models is trained, we can evaluate the model's performance
# 
# first, let's have a look at two generic images downloaded from google. we
# define a image prediction function that does some image preprocessing (so that
# the data fit the data format with which the network was trained) and and
# returns the class with the highest probability (together with the probability
# itself)
pred_img <- function(path) {
  img <- image_load(path, target_size = c(img_width, img_height))
  x <- image_to_array(img)
  x <- x/255
  x <- array_reshape(x, c(1, dim(x)))
  pred <- model %>% predict(x)
  return(c(which(train_generator$class_indices == which.max(pred) - 1), # python: index starts with 0! --> -1
           max(pred)))
}

# also create a helper function to plot an image
img_plot <- function(path){
  img <- image_load(path, target_size = c(img_width, img_height))
  x <- image_to_array(img)
  x <- x/255
  grid::grid.raster(x)
}

img_plot("cocacola.jpg")
pred_img("cocacola.jpg")

img_plot("pepsi.jpg")
pred_img("pepsi.jpg")

img_plot("heineken.jpg")
pred_img("heineken.jpg")


# performance on official test set
test_generator <- flow_images_from_directory(
  test_directory, generator = datagen_val,
  target_size = c(img_width, img_height),
  class_mode = "binary", batch_size = batch_size,
  seed = myseed)

test_performance <- model %>% evaluate_generator(test_generator, steps = 5)
print(paste0("Test accuracy: ", round(test_performance$acc*100,4), "%"))

