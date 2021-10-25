def add_legend(ax, loc=0, fontsize=12, frameon=False):
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
        if "NBm" in l:
            sorts[h, l] = "000" + sort_label
        elif "Universe Machine" in l:
            sorts[h, l] = "zzz" + sort_label
        else:
            sorts[h, l] = sort_label

    # then actually do the sorting
    handles, labels = zip(
        *sorted([hl for hl in zip(handles, labels)], key=lambda hl: sorts[hl])
    )

    ax.legend(
        handles=handles,
        labels=labels,
        loc=loc,
        fontsize=fontsize,
        frameon=frameon,
    )
