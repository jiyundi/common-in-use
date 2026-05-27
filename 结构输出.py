import inspect
import spec_utils  # 你的脚本名

print(f"{spec_utils.__name__}.py")

# 1. 获取并打印所有独立函数 (Functions)
print("├── [独立函数 (Functions)]")
functions = inspect.getmembers(spec_utils, inspect.isfunction)
if functions:
    for name, _ in functions:
        print(f"│    ├── {name}()")
else:
    print("│    └── (无独立函数)")

# 2. 获取所有类 (Classes)
print("├── [类与方法 (Classes & Methods)]")
classes = inspect.getmembers(spec_utils, inspect.isclass)

if classes:
    for class_name, class_obj in classes:
        print(f"│    ├── Class: {class_name}")
        
        # 获取该类内部的所有方法
        # inspect.ismethod 用于类方法，inspect.isfunction 用于类中的普通实例方法
        methods = inspect.getmembers(class_obj, predicate=lambda x: inspect.isfunction(x) or inspect.ismethod(x))
        
        for method_name, _ in methods:
            # 过滤掉 Python 自带的双下划线魔法方法（如 __init__），如果需要看可以删掉这行 filter
            if not method_name.startswith("__"):
                print(f"│    │    ├── {method_name}()")
else:
    print("│    └── (无定义类)")
