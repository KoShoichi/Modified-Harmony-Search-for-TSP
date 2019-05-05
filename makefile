#makefile

CPP = gcc 
INCLUDE_PATH = -I./  

CPPDEFS =-Wall -O2

OBJ = cpu_time.o hstsp.o StdAfx.o main.o

all:main clean
	@echo "all is updated."
	
main:${OBJ}
	${CPP} ${OBJ} -o  $@ -lm 
		
cpu_time.o:cpu_time.c
	${CPP} ${CPPDEFS} -c $< -o $@ ${INCLUDE_PATH}
hstsp.o:hstsp.c
	${CPP} ${CPPDEFS} -c $< -o $@ ${INCLUDE_PATH}
StdAfx.o:StdAfx.c
	${CPP} ${CPPDEFS} -c $< -o $@ ${INCLUDE_PATH}
main.o:main.c
	${CPP} ${CPPDEFS} -c $< -o $@ ${INCLUDE_PATH}
clean:
	rm ${OBJ} 
