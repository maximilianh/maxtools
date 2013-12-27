import os
for subdir, dirs, files in os.walk("."):
    for f in files:
        print os.path.splitext(f)[1]
