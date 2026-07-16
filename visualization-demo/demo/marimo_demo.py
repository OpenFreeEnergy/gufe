"""marimo demo for the gufe framejs visualizations.

Run with:  just dev  (marimo → http://localhost:2718)

marimo is where the py3Dmol pain point shows: py3Dmol does not produce an
anywidget and breaks here, whereas MetaframeWidget (what gufe `.view()` returns)
renders natively.
"""

import marimo

__generated_with = "0.23.13"
app = marimo.App()


@app.cell
def _():
    # The widget the gufe `.view()` integration is built on — renders in marimo
    # with no extra wiring.
    from metaframe_widget import MetaframeWidget

    w = MetaframeWidget.from_code(
        """
        export function onInputs() {
          root.innerHTML = "<div style='font-family:sans-serif;padding:1rem'>"
            + "<h3>framejs in marimo ✓</h3>"
            + "<p>If you can read this, anywidget renders here.</p></div>";
        }
        """,
        height="200px",
    )
    w
    return


@app.cell
def _():
    # The real deliverable: gufe LigandNetwork.view() rendering natively in marimo.
    from pathlib import Path

    import gufe

    data = Path(gufe.__file__).parent / "tests" / "data"
    net = gufe.LigandNetwork.from_graphml((data / "ligand_network.graphml").read_text())
    net.view()
    return


if __name__ == "__main__":
    app.run()
