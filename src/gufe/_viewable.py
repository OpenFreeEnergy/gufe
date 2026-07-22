# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""The ``.view()`` / auto-display mixin shared by every framejs-visualizable object.

Mixing :class:`FramejsViewable` into a gufe class is the *entire* opt-in: the
object then renders as an interactive framejs.io widget in Jupyter/marimo/VSCode
(bare cell or explicit ``.view()``), and :func:`gufe.visualization.framejs.build_cli_url`
can turn it into a shareable browser URL.

Which viz it gets — the frame that draws it and the serializer that feeds that
frame — is looked up by class name in ``framejs.VIZ_REGISTRY``, which walks the
MRO, so a subclass inherits its parent's viz unless it registers its own.

Everything degrades gracefully. With no ``viz`` extra installed, both paths fall
back to the object's ``_legacy_view()`` (RDKit / py3Dmol) if it defines one:
``.view()`` returns it directly, and ``_repr_mimebundle_`` returns its mimebundle
so a bare cell keeps the picture it had before framejs. An object with no legacy
renderer falls through to the plain ``repr``.

Display belongs to this mixin, via ``_repr_mimebundle_`` only. A class must not
also define ``_ipython_display_`` — IPython checks that hook *first* and
short-circuits on it, which would make a bare cell and ``.view()`` disagree.

This module deliberately imports nothing from gufe at import time — the framejs
machinery is pulled in lazily inside each method so ``import gufe`` stays cheap
and dependency-free.
"""

from __future__ import annotations


class FramejsViewable:
    """Mixin giving a gufe object an interactive framejs.io view."""

    def view(self, *, width: str | None = None, height: str | None = None, **opts):
        """Return an interactive framejs.io widget for this object.

        Renders inline in Jupyter, marimo and VSCode (via ``metaframe-widget``).
        Requires the ``viz`` extra (``pip install gufe[viz]``); if it is not
        available this warns and returns the object's ``_legacy_view()``, or
        ``None`` when it has none.

        Parameters
        ----------
        width, height
            CSS sizes for the widget; default to the viz's own preferred size.
        **opts
            Forwarded to the legacy renderer if framejs is unavailable.

        Returns
        -------
        metaframe_widget.MetaframeWidget or None
        """
        from .visualization.framejs import view_object

        return view_object(self, width=width, height=height, **opts)

    def _repr_mimebundle_(self, include=None, exclude=None):
        """Auto-display the interactive framejs widget in notebook front-ends.

        When the ``viz`` extra is not installed or framejs is otherwise
        unavailable, returns the legacy renderer's mimebundle if the object has
        one, else ``None`` so the notebook prints the plain ``repr``. Silent
        either way — this runs on every display of the object.
        """
        from .visualization.framejs import repr_mimebundle

        return repr_mimebundle(self, include=include, exclude=exclude)
