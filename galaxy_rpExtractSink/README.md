# Galaxy rpExtractSink

Galaxy tool that takes for input an SBML file and a compartment ID and extracts the chemical species and formats then to a CSV valid as an input sink for RetroPath2.0

## Getting Started

This is a docker galaxy tools, and thus, the docker needs to be built locally where Galaxy is installed. 

### Prerequisites

TODO

### Installing

Create a new section in the Galaxy tool_cong.xml from the config file:

```
<section id="retro" name="Retro Nodes">
  <tool file="/local/path/config_rpExtractSink.xml" />
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
docker build -t ibisba/rpextractsink .
```

Make sure that the galaxy instance contains the following in its job_conf.xml:

```
<destination id="docker_data_readonly" runner="local">
	<param id="docker_enabled">true</param>
	<param id="docker_sudo">false</param>
	<param id="docker_auto_rm">true</param>
	<param id="docker_set_user">root</param>
  <param id="docker_run_extra_arguments">--mount source=rpcache,destination=/home/rpInspect/cache,readonly --mount source=component_contribution_data,destination=/home/rpInspect/component_contribution/data,readonly</param>

</destination>
```

and 

```
<tool id="rpExtractSink" destination="docker_data_readonly" />
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
