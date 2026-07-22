# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""The ``.view()`` / auto-display mixin shared by every framejs-visualizable object.

Mixing :class:`FramejsViewable` into a gufe class is the *entire* opt-in: the
object then renders as an interactive framejs.io widget in Jupyter/marimo/VSCode
(bare cell or explicit ``.view()``), and :func:`gufe.visualization.framejs.build_cli_url`
can turn it into a shareable browser URL.

Which viz and which serialization it gets are looked up by class name in
``framejs.CANONICAL_VIZ`` / ``framejs._PAYLOAD_BUILDERS`` — both walk the MRO, so a
subclass inherits its parent's viz unless it registers its own.

Everything degrades gracefully: with no ``viz`` extra installed, ``.view()`` falls
back to the object's ``_legacy_view`` (RDKit / py3Dmol) if it defines one, and
``_repr_mimebundle_`` returns ``None`` so the notebook prints the plain ``repr``.

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
        available this returns the object's legacy view, or ``None`` when it has
        none.

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

        Returns ``None`` (so the notebook falls back to the plain ``repr``) when
        the ``viz`` extra is not installed or framejs is otherwise unavailable.
        """
        from .visualization.framejs import repr_mimebundle

        return repr_mimebundle(self, include=include, exclude=exclude)
