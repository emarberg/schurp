

def cached_value(cache):
    def wrap(f):
        @classmethod
        def wrapped_f(cls, *args):
            key = tuple(args)
            if key not in cache:
                cache[key] = f(cls, *args)
            return cache[key]
        return wrapped_f
    return wrap
