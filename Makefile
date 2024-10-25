genAddr:
	g++ -Wno-write-strings genAddrBTC.cpp -o genAddr
clean:
	@rm -f genAddr
	@echo clean done !