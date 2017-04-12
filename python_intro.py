def sphere_volume(r):
    volume = (4./3)*3.14159 * r**3
    return volume

def first_half(r):
    length = len(r)
    half = length / 2
    int_half = int(half)
    a_kisss = r[0:int_half]
    return a_kisss

def backward(s):
    leng = len(s)
    a_kiss = s[leng::-1]
    return a_kiss

def list_ops():
    my_list = ["bear","ant","dog","cat"]
    my_list.append("eagle")
    my_list[2] = "fox"
    my_list.remove(my_list[1])
    my_list.sort(key=str.lower)
    my_list = backward(my_list)
    return my_list

def pig_latin(r):
    vowel = ['a', 'e', 'i', 'o', 'u']
    fir = r[0]
    if fir in vowel:
        rod = r + "hay"
    else:
        new = r[1:]
        rod = new + fir + "ay"
    return rod


def palindrome():
	b = 0
	c = 100
	d = 100
	while c < 1000:
		while d < 1000:
			a = c * d
			aa = str(a)
			leng = len(aa)
			zz = aa[leng::-1]
			z = int(zz)
			d = d + 1
			if a == z and a > b:
				b = a
		c = c + 1
		d = 100
	return b

def alt_harmonic(i):
	series = [((-1)**(n+1))/float(n) for n in range(1,i+1)]
	result = sum(series)
	return result

if __name__ == "__main__":
    print("Hello, world!")
    print sphere_volume(5)
