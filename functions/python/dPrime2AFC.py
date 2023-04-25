def dprime_2afc(correct_side, response_side):
  """Calculate the d' score for a 2AFC task.

  Args:
    correct_side: The assigned correct side.
    response_side: The side chosen by the animal (NaN tolerated).

  Returns:
    The d' score, currently without 1/2^(1/2) correction.
  """

  # Input check
  correct_side = correct_side[~np.isnan(response_side)]
  response_side = response_side[~np.isnan(response_side)]

  # Determine the side to look at
  side_codes = np.unique(response_side)
  the_side = side_codes[0]

  hit_rate = np.sum(correct_side == the_side & response_side == the_side) / np.sum(correct_side == the_side)
  false_alarm_rate = np.sum(correct_side != the_side & response_side == the_side) / np.sum(correct_side != the_side)

  d_prime = scipy.stats.norm.ppf(hit_rate) - scipy.stats.norm.ppf(false_alarm_rate)

  # Handle degenerate cases
  if d_prime == np.inf or d_prime == -np.inf:
    if hit_rate != 1 and false_alarm_rate != 0:
      if hit_rate == 1:
        d_prime = scipy.stats.norm.ppf(1 - false_alarm_rate)
      elif false_alarm_rate == 1:
        d_prime = scipy.stats.norm.ppf(1 - hit_rate)
      else:
        d_prime = 0
    elif hit_rate - false_alarm_rate < 0.5:
      d_prime = 0

  return d_prime
