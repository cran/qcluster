# Predict Clustering using gmix Model
#
# This function predicts the clustering for new data based on a fitted model of class `mbcfit`.
#
# Parameters:
# object - An object of class `mbcfit` representing the fitted model.
# newdata - A data frame or matrix containing new data for prediction.
#
# Returns:
# A vector indicating the cluster assignment for each observation in the new data.
#
# Example:
# Assuming `model` is a fitted object of class `mbcfit`
# and `new_data` is the data to be clustered
# predicted_clusters <- predict.gmix(model, new_data)
# print(predicted_clusters)

predict.mbcfit <- function(object, newdata, ...) {
  xname <- deparse(substitute(object))

  if (!inherits(object, "mbcfit")) {
    stop(xname, ' is not of class "mbcfit"')
  } else if (!object$info$code %in% c(1, 2)) {
    print(object)
    stop(xname, " does not contain enough information to estimate new clustering")
  }

  # Check data input
  newdata <- .ckdat(newdata)

  if (newdata$p != object$P) {
    stop(sprintf("'data' dimensionality (%d) is non conformable with that of '%s' (%d)", newdata$p, xname, object$P))
  }

  cluster <- .dpclassify(newdata$data, object$params)
  return(cluster)
}
