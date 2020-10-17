SUBDIRS = src

distdir = $(PACKAGE)

dist-hook:
	find $(distdir) -type f -name '*~' -exec rm {} \;

install:
	mv src/eprem ../../bin/eprem
