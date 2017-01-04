import pyautogui as gui
import time
import sys

def main():
	input1 = sys.argv[1]
	print_to_kodi(input1)

def print_to_kodi(s):
	nt = 3
	for i in range(nt):
		ss = 'time before writing in kodi is ' + str(nt-i) + 's'
		print(ss)
		time.sleep(1)
	gui.typewrite(s, interval=0.25)


if __name__ == "__main__":
	main()
