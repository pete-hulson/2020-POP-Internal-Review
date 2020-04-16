ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: pop 

ms: ms.tpl 
	@admb -f ms.tpl 

pop: pop.tpl 
	@admb  pop.tpl 

clean:
	@$(RM) *.cpp
	@$(RM) *.htp
	@$(RM) *.obj
