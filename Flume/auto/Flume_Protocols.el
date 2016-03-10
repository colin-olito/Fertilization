(TeX-add-style-hook
 "Flume_Protocols"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "amsmath"
    "amsfonts"
    "amssymb"
    "inputenc"
    "color"
    "soul"
    "times"
    "fullpage")))

