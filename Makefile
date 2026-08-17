# CaRaSh local spherical shell experiment (Microsoft NMake)

CC = cl
CFLAGS = /nologo /std:c11 /O2 /W4 /D "NOMINMAX"

!IFNDEF L
L = 17
!ENDIF

TARGET = carash.exe

all: $(TARGET)

$(TARGET): carash.c
	$(CC) $(CFLAGS) /D "L=$(L)" carash.c /Fe:$(TARGET)

run: $(TARGET)
	$(TARGET)

clean:
	-del $(TARGET) *.obj
