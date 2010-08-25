#!/bin/bash

# FIXME: Do we really need this file?
#
# Here's a better fix:
#
#;; store backup files in /tmp
#(setq backup-directory-alist
#      `((".*" . ,temporary-file-directory)))
#(setq auto-save-file-name-transforms
#      `((".*" ,temporary-file-directory t)))

find . -name '*~' | xargs rm -f
