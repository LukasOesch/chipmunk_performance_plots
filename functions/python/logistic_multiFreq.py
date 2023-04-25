def logistic_multifreq(performance, binary_response_side, stim_freq, category_boundary):
  """
  Calculates the logistic model for multiple frequencies.

  Args:
    performance: A vector of performance values.
    binary_response_side: A vector of binary response values.
    stim_freq: A vector of stimulus frequencies.
    category_boundary: The category boundary.

  Returns:
    A tuple of four values:
      * The choice model statistics.
      * The x-data for the psychometric function.
      * The y-data for the psychometric function.
      * The frequency bias.
      * The frequency sensitivity.
  """

  # Create the choice predictors.
  choice_predictors = []
  side_chosen = []
  for n in range(len(performance)):
    if not np.isnan(performance[n]):
      choice_predictors.append(stim_freq[n] - category_boundary)
      side_chosen.append(binary_response_side[n])

  # Fit the logistic model.
  choice_coef, _, choice_model_stats = statsmodels.api.glm.GLM.fit(
      choice_predictors, side_chosen, family=statsmodels.api.families.Binomial(), link='logit')

  # Calculate the x-data for the psychometric function.
  x_data = np.arange(min(stim_freq) - category_boundary, max(stim_freq) - category_boundary, 0.5)

  # Calculate the y-data for the psychometric function.
  p_data = 1. / (1 + np.exp(-(choice_coef[0] + choice_coef[1] * x_data)))

  # Calculate the frequency bias.
  freq_bias = category_boundary + (np.log(1) - choice_coef[0]) / choice_coef[1]

  # Calculate the frequency sensitivity.
  freq_sensitivity = abs(choice_coef[1]) / 4

  return choice_model_stats, x_data, p_data, freq_bias, freq_sensitivity
