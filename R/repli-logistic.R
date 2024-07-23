# Logistic Distribution - on tanh scale - link functions.
# These are the Logistic Distribution CDF. The tanh/arctanh functions would
# actually rescale the predicted values before link function as well as
# afterwards. However, the range of our responses will actually follow the
# range of tanh - scale Logistic probabilities to be on the range (-1,1), so
# that they look good in IGV.
qlogistanh <- \(q) qlogis(0.5*q + 0.5)
plogistanh <- \(x) (2*plogis(x) - 1)
