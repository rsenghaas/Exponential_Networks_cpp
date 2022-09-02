numbers = [30,16,4,2,2,5,25,1,5,5]


for i in range(2**len(numbers)):
    signs = []
    temp_i = i
    for _ in range(len(numbers)):
        if temp_i % 2 == 0:
            signs.append(-1);
        else:
            signs.append(1);
        temp_i -= (temp_i % 2);
        temp_i /= 2;
    total = 0
    for k in range(len(numbers)):
        total += signs[k] * numbers[k]

    if abs(total - 22) <= 0:
        sret = ""
        for k in range(len(numbers)):
            if signs[k] < 0:
                sret += f"-{numbers[k]}"
            elif k != 0:
                sret += f"+{numbers[k]}"
            else:
                 sret += f"{numbers[k]}"
        sret += f"={total}"
        print(sret)


            
