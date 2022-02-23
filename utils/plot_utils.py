from pathlib import Path

import numpy as np
from matplotlib import ticker


# ======================================================================================
#
# functions for commonly used plotting ideas
#
# ======================================================================================
def get_sim_dirname(sim_loc):
    sim_path = Path(sim_loc)
    if sim_path.name == "run":
        sim_path = sim_path.parent

    # now I should have the directory
    return sim_path.name


def plot_label(sim, sim_share_type, axis_name, include_z=True):
    label = sim.names[axis_name]
    # and include the redshift if it's different for each sim
    if sim_share_type == "last" and include_z:
        z = 1 / sim.ds.scale_factor - 1
        # don't include if at z=1.5
        if not 1.49 < z < 1.51:
            label += f": z = {z:.1f}"

    return label


def add_legend(ax, loc=0, fontsize=12, frameon=False, **kwargs):
    # make a default legend
    ax.legend()
    # then get everything from it
    handles, labels = ax.get_legend_handles_labels()

    if len(handles) == 0:
        return

    # then figure out how to sort them
    sorts = dict()
    for h, l in zip(handles, labels):
        # get rid of the redshift part of the label
        if ": z = " in l:
            sort_label = l[: l.find(":")]
        else:
            sort_label = l
        # add extras to manually change the sort for labels I want at the front
        # or back
        if "NBm" in l or "$f_{boost}=1$, $f_{HN,0}=50$%" in l:
            # ! is the first real ASCII character
            sorts[h, l] = "!!!" + sort_label
        elif "Universe Machine" in l or l.startswith("M < ") or l.startswith("M > "):
            sorts[h, l] = "zzz" + sort_label
        else:
            sorts[h, l] = sort_label

    # then actually do the sorting
    handles, labels = zip(
        *sorted([hl for hl in zip(handles, labels)], key=lambda hl: sorts[hl])
    )

    return ax.legend(
        handles=handles,
        labels=labels,
        loc=loc,
        fontsize=fontsize,
        frameon=frameon,
        **kwargs,
    )


def moving_percentiles(xs, ys, percentile, dx, bins=None):
    if bins is None:
        bins = np.arange(min(xs), max(xs) + dx, dx)

    bin_centers = []
    radii_percentiles = []
    for idx in range(len(bins) - 1):
        lower = bins[idx]
        upper = bins[idx + 1]

        # then find all clusters in this mass range
        mask_above = xs > lower
        mask_below = xs < upper
        mask_good = np.logical_and(mask_above, mask_below)

        good_ys = ys[mask_good]
        if len(good_ys) > 1:
            radii_percentiles.append(np.percentile(good_ys, percentile))
            # the bin centers will be the mean in log space
            bin_centers.append(np.mean([lower, upper]))

    return np.array(bin_centers), np.array(radii_percentiles)


def shaded_region(
    ax,
    xs,
    ys,
    color,
    p_lo=10,
    p_hi=90,
    dx=0.2,
    alpha=0.5,
    zorder=0,
    cut_x=np.inf,
    log_x=False,
    label=None,
    bins=None,
):
    # cut_x is the x value above which we plot dashed lines.
    if log_x:
        xs = np.log10(xs)
        # don't take the log of cut, since it's not used until later.

    c_lo, p_lo = moving_percentiles(xs, ys, p_lo, dx, bins)
    c_50, p_50 = moving_percentiles(xs, ys, 50, dx, bins)
    c_hi, p_hi = moving_percentiles(xs, ys, p_hi, dx, bins)

    if log_x:
        c_lo = 10 ** c_lo
        c_50 = 10 ** c_50
        c_hi = 10 ** c_hi

    # The X values should be the same for all percentiles
    assert np.allclose(c_lo, c_50)
    assert np.allclose(c_hi, c_50)

    plot_line_with_cut(
        ax, c_50, p_50, cut_x, color=color, zorder=zorder + 1, label=label
    )
    ax.fill_between(x=c_lo, y1=p_lo, y2=p_hi, color=color, alpha=alpha, zorder=zorder)


def plot_line_with_cut(ax, xs, ys, cut_x, label, ls1="-", ls2="--", **kwargs):
    """
    plot a line, with different linestyles before and after a given point
    """
    # find the range at which we need to start plotting dashed lines
    # the below code breaks if the cut is below all values, so do that independently
    if cut_x < np.min(xs):
        ax.plot(xs, ys, ls=ls2, label=label, **kwargs)
        return
    # somewhere in the middle
    for cut_idx in range(len(xs)):
        if xs[cut_idx] > cut_x:
            # cut at the index we're at now
            ax.plot(xs[:cut_idx], ys[:cut_idx], ls=ls1, label=label, **kwargs)
            # then the unreliable range. We start at the place the unreliable
            # range ended to get a continuous line
            ax.plot(xs[cut_idx - 1 :], ys[cut_idx - 1 :], ls=ls2, label=None, **kwargs)
            break  # end the loop so we don't plot more
    else:  # no break found
        ax.plot(xs, ys, ls=ls1, label=label, **kwargs)


# Function to use to set the ticks
@ticker.FuncFormatter
def nice_log_formatter(x, pos):
    exp = np.log10(x)
    # this only works for labels that are factors of 10. Other values will produce
    # misleading results, so check this assumption.
    assert np.isclose(exp, int(exp))

    # for values between 0.01 and 100, just use that value.
    # Otherwise use the log.
    if abs(exp) < 3:
        return f"{x:g}"
    else:
        return "$10^{" + f"{exp:.0f}" + "}$"


def nice_log_axis(ax, which):
    if which not in ["x", "y", "both"]:
        raise ValueError("error specifying axis in nice_log_axis")
    if which in ["x", "both"]:
        ax.xaxis.set_major_formatter(nice_log_formatter)
    if which in ["y", "both"]:
        ax.yaxis.set_major_formatter(nice_log_formatter)


# ======================================================================================
#
# functions for KDE histograms
#
# ======================================================================================
def gaussian(x, mean, variance):
    """
    Normalized Gaussian Function at a given value.

    Is normalized to integrate to 1.

    :param x: value to calculate the Gaussian at
    :param mean: mean value of the Gaussian
    :param variance: Variance of the Gaussian distribution
    :return: log of the likelihood at x
    """
    exp_term = np.exp(-((x - mean) ** 2) / (2 * variance))
    normalization = 1.0 / np.sqrt(2 * np.pi * variance)
    return exp_term * normalization


def kde(x_grid, x_values, width, weights, log=False):
    ys = np.zeros(x_grid.size)

    if log:
        x_grid = np.log10(x_grid)
        x_values = np.log10(x_values)

    for x, weight in zip(x_values, weights):
        ys += weight * gaussian(x_grid, x, width ** 2)

    # normalize the y value
    ys = np.array(ys)
    ys = 10 * ys / np.sum(ys)
    return ys
