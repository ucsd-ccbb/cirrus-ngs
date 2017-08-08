def test_gen():
    for i in range(5):
        yield i, i+1

result = ["proj"]
for x in test_gen():
    print(result + list(x))

print(result)
