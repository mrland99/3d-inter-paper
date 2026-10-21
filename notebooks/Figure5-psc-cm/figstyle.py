"""Shared matplotlib style for manuscript figure panels: Arial + Illustrator-editable text.

Usage (in a panel script):
    import figstyle; figstyle.apply()      # call once, before plotting
"""
import matplotlib, matplotlib.font_manager as fm


def apply(base=11):
    """Arial everywhere + vector text that stays editable in Illustrator."""
    assert "Arial" in {f.name for f in fm.fontManager.ttflist}, \
        "Arial not found by matplotlib (install arial.ttf to ~/.local/share/fonts and clear the mpl cache)"
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        # make inline math ($...$) Arial too, so it matches the plain text
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
        # keep text as text (NOT outlines) so Illustrator can edit it:
        "svg.fonttype": "none",     # SVG -> <text> elements
        "pdf.fonttype": 42, "ps.fonttype": 42,   # embed TrueType, selectable/editable
        "font.size": base, "axes.titlesize": base + 1, "axes.labelsize": base,
        "xtick.labelsize": base - 1, "ytick.labelsize": base - 1, "legend.fontsize": base - 1,
    })
