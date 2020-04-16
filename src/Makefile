ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: pop 

pop: pop.tpl 
	@admb  pop.tpl 
ifdef ComSpec
	$(COPY) pop.exe mypath
else
	$(COPY) pop ~/
endif

clean:
	@$(RM) *.cpp
	@$(RM) *.htp
	@$(RM) *.obj
