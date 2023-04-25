def calculate_wilson_score_interval(data):
  """Calculate the Wilson score confidence interval for binomial distributions.

  This confidence interval is asymmetric, bounded by 0 and 1, and "pulls" towards 0.5.
  Currently, this function only handles an alpha level of 0.05.

  Args:
    data: A vector of 0s and 1s containing the observed outcomes.

  Returns:
    A tuple of two floats, the lower and upper bounds of the confidence interval.
  """

  # Take the average of all the specified values.
  p = np.nanmean(data)

  # Count the number of valid observations.
  n = np.sum(~np.isnan(data))

  # The critical value for an alpha level of 0.05.
  z = 1.95996

  # Calculate the lower bound of the confidence interval.
  lowerBound = (p + z**2 / (2 * n) - z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / (1 + z**2 / n)

  # Calculate the upper bound of the confidence interval.
  upperBound = (p + z**2 / (2 * n) + z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / (1 + z**2 / n)

  return lowerBound, upperBound
