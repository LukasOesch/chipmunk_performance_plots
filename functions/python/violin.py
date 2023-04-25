import numpy as np
import matplotlib.pyplot as plt

class Violin(object):
    """
    Creates violin plots for some data

    Args:
        data: A data vector.
        pos: The X-axis position of the violin plot.
        width: The width of the violin plot in axis space.
        bandwidth: The bandwidth of the kernel density estimation.
        violin_color: The fill color of the violin area and data points.
        violin_alpha: The transparency of the violin area and data points.
        edge_color: The color of the violin area outline.
        box_color: The color of the box, whiskers, and the outlines of the median point and the notch indicators.
        box_width: The width of the box between the quartiles in axis space (default 10% of Violin plot width, 0.03)
        median_color: The fill color of the median and notch indicators.
        show_data: Whether to show data points.
        show_notches: Whether to show notch indicators.
        show_mean: Whether to show mean indicator.

    """

    def __init__(self, data, pos, width, bandwidth, violin_color, violin_alpha, edge_color, box_color, box_width, median_color, show_data, show_notches, show_mean):
        self.data = data
        self.pos = pos
        self.width = width
        self.bandwidth = bandwidth
        self.violin_color = violin_color
        self.violin_alpha = violin_alpha
        self.edge_color = edge_color
        self.box_color = box_color
        self.box_width = box_width
        self.median_color = median_color
        self.show_data = show_data
        self.show_notches = show_notches
        self.show_mean = show_mean

        # Calculate the mean and standard deviation of the data.
        self.mean = np.mean(self.data)
        self.std = np.std(self.data)

        # Calculate the kernel density estimate of the data.
        self.kde = np.lib.stats.gaussian_kde(self.data)

        # Plot the violin.
        self.violin_plot = plt.fill([self.pos + width / 2 * self.kde(x) for x in np.linspace(min(self.data), max(self.data))], [x for x in np.linspace(min(self.data), max(self.data))], self.violin_color, alpha=self.violin_alpha)
        self.violin_plot.set_edgecolor(self.edge_color)

        # Plot the box.
        self.box_plot = plt.fill([self.pos - self.box_width / 2, self.pos + self.box_width / 2], [self.quantile(self.data, 0.25), self.quantile(self.data, 0.75)], self.box_color, alpha=self.violin_alpha)
        self.box_plot.set_edgecolor(self.edge_color)

        # Plot the whiskers.
        self.whisker_plot = plt.plot([self.pos - self.box_width / 2, self.pos - self.box_width / 2], [self.quantile(self.data, 0.25) - 1.5 * self.std, self.quantile(self.data, 0.25)], self.box_color, alpha=self.violin_alpha)
        self.whisker_plot = plt.plot([self.pos + self.box_width / 2, self.pos + self.box_width / 2], [self.quantile(self.data, 0.75) + 1.5 * self.std, self.quantile(self.data, 0.75)], self.box_color, alpha=self.violin_alpha)

        # Plot the median.
        self.median_plot = plt.plot([self.pos, self.pos], [self.quantile(self.data, 0.5), self.quantile(self.data, 0.5
        self.median_plot.set_color(self.median_color)
        self.median_plot.set_linewidth(2)

        # Plot the mean.
        if self.show_mean:
            self.mean_plot = plt.plot([self.pos, self.pos], [self.mean, self.mean], self.median_color, linewidth=2)

        # Plot the data points.
        if self.show_data:
            self.data_plot = plt.scatter(self.pos + np.random.uniform(-0.5, 0.5, size=len(self.data)), self.data, s=50, alpha=self.violin_alpha)

        # Set the axis limits.
        plt.xlim([self.pos - self.width / 2, self.pos + self.width / 2])
        plt.ylim([min(self.data), max(self.data)])

        # Show the plot.
        plt.show()

    def set_edge_color(self, color):
        self.violin_plot.set_edgecolor(color)

    def get_edge_color(self):
        return self.violin_plot.get_edgecolor()

    def set_median_color(self, color):
        self.median_plot.set_color(color)

    def get_median_color(self):
        return self.median_plot.get_color()

