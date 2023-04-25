def plot_areaerrorbar(data, options):
  """
  Plots the mean and standard deviation of a set of data filling the space
  between the positive and negative mean error using a semi-transparent
  background, completely customizable.

  Args:
    data: A data matrix, with rows corresponding to observations and
      columns to samples.
    options: An optional struct that contains the customized parameters.
      * options.handle: The figure handle to plot the result.
      * options.color_area: The RGB color of the filled area.
      * options.color_line: The RGB color of the mean line.
      * options.alpha: The alpha value for transparency.
      * options.line_width: The mean line width.
      * options.x_axis: The X time vector.
      * options.error: The type of error to plot (+/-).
        * If 'std', one standard deviation.
        * If 'sem', the standard error of the mean.
        * If 'var', one variance.
        * If 'c95', a 95% confidence interval.

  Returns:
    None.
  """

  # Default options
  if len(options) == 0:
    options = {}
    options['handle'] = plt.figure(1)
    options['color_area'] = [128, 193, 219] / 255  # Blue theme
    options['color_line'] = [52, 148, 186] / 255
    options['alpha'] = 0.5
    options['line_width'] = 2
    options['error'] = 'std'

  # Check if the x-axis is specified.
  if 'x_axis' not in options:
    options['x_axis'] = np.arange(1, len(data) + 1)

  # Compute the mean and standard deviation of the data matrix.
  data_mean = np.nanmean(data, axis=1)
  data_std = np.nanstd(data, axis=0)

  # Type of error plot.
  switch(options['error']):
    case 'std':
      error = data_std
    case 'sem':
      error = data_std / np.sqrt(data.shape[1])
    case 'var':
      error = data_std ** 2
    case 'c95':
      error = data_std / np.sqrt(data.shape[1]) * 1.96

  # Plot the result.
  x_vector = np.concatenate((options['x_axis'], np.fliplr(options['x_axis'])))
  patch = plt.fill(x_vector, data_mean + error, options['color_area'])
  patch.set_edgecolor('none')
  patch.set_facealpha(options['alpha'])
  plt.hold(True)
  plt.plot(options['x_axis'], data_mean, 'color', options['color_line'], ...
           'linewidth', options['line_width'])
  plt.hold(False)

  return None
