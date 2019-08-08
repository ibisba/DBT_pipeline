# Galaxy rpMergeSBML

Galaxy tool that takes a tar.xz with a list of rpSBML models and merges them with another SBML file

## Getting Started

This is a docker galaxy tools, and thus, the docker needs to be built locally where Galaxy is installed. 

### Prerequisites

TODO

### Installing

Create a new section in the Galaxy tool_cong.xml from the config file:

```
<section id="retro" name="Retro Nodes">
  <tool file="/local/path/galaxyRPmerge.xml" />
</section>
```

Make sure that docker can be run in root:

```
sudo groupadd docker
sudo gpasswd -a $USER docker
sudo service docker restart
```

Build the docker image:

```
docker build -t ibisba/rpmergesbml .
```

Make sure that the following job_conf.xml looks like this:

```
<?xml version="1.0"?>
<job_conf>
  <plugins>
    <plugin id="local" type="runner" load="galaxy.jobs.runners.local:LocalJobRunner" workers="4"/>
  </plugins>
  <destinations default="docker_local">
    <destination id="local" runner="local" />
    <destination id="docker_local" runner="local">
      <param id="docker_enabled">true</param>
      <param id="docker_sudo">false</param>
      <param id="docker_auto_rm">true</param>
      <param id="docker_set_user">root</param>
    </destination>
  </destinations>
</job_conf>
```

It is important to run the docker as root user since we will be calling a script that writes files to a temporary folder inside the docker before sending bask to Galaxy

## Running the tests

TODO

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Galaxy](https://galaxyproject.org) - The Galaxy project

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

TODO

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson