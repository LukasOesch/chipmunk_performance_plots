def violinplot(data, cats, **kwargs):
    """
    Violinplots plots violin plots of some data and categories

    Args:
        data: A data vector or a 2D data matrix.
        cats: A vector of categories, or None if there are no categories.
        **kwargs: Keyword arguments to be passed to the `Violin` constructor.

    Returns:
        An array of `Violin` objects.
    """

    # Check the input arguments.
    if not isinstance(data, (list, np.ndarray)):
        raise TypeError('data must be a list or a NumPy array')
    if cats is not None and not isinstance(cats, (list, np.ndarray)):
        raise TypeError('cats must be a list or a NumPy array')

    # If data is a 2D data matrix, convert it to a list of 1D data vectors.
    if len(data.shape) == 2:
        data = data.T

    # Create an array of `Violin` objects.
    violins = []
    for i, d in enumerate(data):
        if cats is None:
            violins.append(Violin(d, i, **kwargs))
        else:
            violins.append(Violin(d, cats[i], **kwargs))

    return violins

