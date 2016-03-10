(TeX-add-style-hook
 "Olito_Marshall_Connallon_AmNat_two_sex_spawning_time_R0"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("lineno" "left")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "lineno"
    "titlesec"
    "amsmath"
    "amsfonts"
    "amssymb"
    "color"
    "soul"
    "times"
    "fullpage")))

