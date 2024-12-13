all:
	g++ -Wno-write-strings genAddrBTC.cpp -o hiiu
clean:
	@rm -f hiiu
	@echo clean done !