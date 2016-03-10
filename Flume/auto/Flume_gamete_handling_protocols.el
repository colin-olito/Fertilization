(TeX-add-style-hook
 "Flume_gamete_handling_protocols"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "titlesec"
    "amsmath"
    "amsfonts"
    "amssymb"
    "color"
    "soul"
    "times"
    "fullpage")))

