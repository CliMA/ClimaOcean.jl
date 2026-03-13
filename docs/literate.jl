using Literate

script_path = ARGS[1]
literated_dir = ARGS[2]

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

@time basename(script_path) Literate.markdown(script_path, literated_dir;
                                              flavor = Literate.DocumenterFlavor(),
                                              execute = true)
