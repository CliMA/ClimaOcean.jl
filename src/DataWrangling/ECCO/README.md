# Setting `ECCO_USERNAME` and `ECCO_PASSWORD` environment variables for downloading ECCO datasets

The first step is to find the username and password for your "programmatic API" credentials on NASA's Earthdrive.
For this you have to either login or make an account via the "EARTHDATA login":

> https://urs.earthdata.nasa.gov

Either register and then sign in or, if you are already registered, sign in. Next, navigate to the ECCO drive:

> https://ecco.jpl.nasa.gov/drive/

This should produce a screen similar to the following:

![image](https://github.com/user-attachments/assets/490d9098-aece-4e9c-82d7-3ec86e833347)

showing your Programmatic API credentials -- except in place of the black boxes that say `your_username` and `cRaZYpASSwORD`,
you should see _your_ username and password.
Copy the content of `Username:` to the environment variable `ECCO_USERNAME` and the content of `Password` to `ECCO_PASSWORD`,
either in a file:

```bash
export ECCO_USERNAME=your_username
export ECCO_PASSWORD=cRaZYpASSwORD
```

or within Julia by writing

```julia
ENV["ECCO_USERNAME"] = "your_username"
ENV["ECCO_PASSWORD"] = "cRaZYpASSwORD"
```
