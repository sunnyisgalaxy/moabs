DIRS=bsmap ExactNumCI mcall mcomp

all:
	@for i in $(DIRS); do \
		make -C $${i}  SRC_ROOT=$(MOABS_ROOT) OPT=1; \
	done;

install:
	@for i in $(DIRS); do \
		make -C $${i}  SRC_ROOT=$(MOABS_ROOT) OPT=1 install; \
	done;

test:
	@for i in $(DIRS); do \
		make -C $${i}  SRC_ROOT=$(MOABS_ROOT) test; \
	done;

clean:
	@for i in $(DIRS); do \
		make -C $${i}  SRC_ROOT=$(MOABS_ROOT) clean; \
	done;
