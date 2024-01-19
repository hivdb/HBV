import matplotlib.pyplot as plt
import numpy as np # noqa


def plot_hist(
        save_path,
        x, y,
        title, x_label, y_label,
        x_color_cut=5):

    f, ax = plt.subplots()

    colors = []
    for v in x:
        if v > x_color_cut:
            colors.append('gray')
        else:
            colors.append('white')

    ax.bar(
        x, y, width=0.7,
        align='center', color=colors, edgecolor='black')

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    ax.ticklabel_format(style='plain')

    ax.set_ylim(ymin=0.1, ymax=max(y) * 1.1)

    ax.get_yaxis().get_major_locator().set_params(integer=True)
    ax.get_xaxis().get_major_locator().set_params(integer=True)

    f.savefig(save_path)


