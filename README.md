<img src="http://i.imgur.com/mGFMeEC.png" alt="It's going to be alright." />
doctor-fishman
================
It's going to be all right.

### Installation
tbd

### Development and Testing
tbd

#### Option 1
As you make changes to your dashboard app, you can run your dashboard app from within RStudio Server using the `shiny` R package:
```
library(shiny)
runApp('path/to/my/app')
```
#### Option 2
Shiny Server also supports per-user applications, which allow you to navigate directly to your app from your browser. Shiny Server expects that user-level apps are stored in `~/ShinyApps`. Assuming you cloned this repo into the root level of your home directory, the following should work:
```
ln -s shiny-dashboards/app ShinyApps
```
You can then navigate to your dashboard:
tbd

### Releasing
tbd

Current release process:
1. Once you're happy with the changes on your local copy, submit a pull request and tag someone on data science who's familiar with Shiny, or even better, the dashboard you've worked on.
2. Get the thumbs up to merge the pull request.
3. ???
4. Profit
