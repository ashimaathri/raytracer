def hexToReal(hexString):
    red = int(hexString[:2], 16) / 255.0
    green = int(hexString[2:4], 16) / 255.0
    blue = int(hexString[4:6], 16) / 255.0
    return (red, green, blue)

print(hexToReal("F1D4AF"))
