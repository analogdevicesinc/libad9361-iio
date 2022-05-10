#!/bin/sh -e

export TRAVIS_API_URL="https://api.travis-ci.com"
LOCAL_BUILD_DIR=${LOCAL_BUILD_DIR:-build}

HOMEBREW_NO_INSTALL_CLEANUP=1
export HOMEBREW_NO_INSTALL_CLEANUP

LIBIIO_BRANCH=master

PYTHON=python3
export PYTHON

# This needs to be duplicated inside 'inside_docker.sh'
# It's the common convention between host & container
INSIDE_DOCKER_BUILD_DIR=/docker_build_dir

# Add here all the common env-vars that should be propagated
# to the docker image, simply by referencing the env-var name.
# The values will be evaluated.
#
# Make sure to not pass certain stuff that are specific to the host
# and not specific to inside-the-docker (like TRAVIS_BUILD_DIR)
#
# If these nothing should be passed, then clear or
#'unset INSIDE_DOCKER_TRAVIS_CI_ENV' after this script is included
INSIDE_DOCKER_TRAVIS_CI_ENV="TRAVIS TRAVIS_COMMIT TRAVIS_PULL_REQUEST OS_TYPE OS_VERSION ARTIFACTNAME PACKAGE_TO_INSTALL"

echo_red()   { printf "\033[1;31m$*\033[m\n"; }
echo_green() { printf "\033[1;32m$*\033[m\n"; }
echo_blue()  { printf "\033[1;34m$*\033[m\n"; }

get_script_path() {
	local script="$1"

	[ -n "$script" ] || return 1

	if [ -f "CI/travis/$script" ] ; then
		echo "CI/travis/$script"
	elif [ -f "ci/travis/$script" ] ; then
		echo "ci/travis/$script"
	elif [ -f "${LOCAL_BUILD_DIR}/$script" ] ; then
		echo "${LOCAL_BUILD_DIR}/$script"
	else
		return 1
	fi
}

get_ldist() {
	case "$(uname)" in
	Linux*)
		. /etc/os-release
		if ! command_exists dpkg ; then
			echo $ID-$VERSION_ID-$(uname -m)
		else
			echo $ID-$VERSION_ID-$(dpkg --print-architecture)
		fi
		;;
	Darwin*)
		echo "darwin-$(sw_vers -productVersion)"
		;;
	*)
		echo "$(uname)-unknown"
		;;
	esac
	return 0
}

__brew_install_if_not_exists() {
	brew ls --versions "$1" || \
		brew install "$1"
}

brew_install_if_not_exists() {
	while [ -n "$1" ] ; do
		__brew_install_if_not_exists "$1" || return 1
		shift
	done
}

prepare_docker_image() {
	local DOCKER_IMAGE="${OS_TYPE}:${OS_VERSION}"
	# If arch is specified, setup multiarch support
	if [ -n "$OS_ARCH" ] ; then
		sudo apt-get -qq update
		sudo DEBIAN_FRONTEND=noninteractive apt-get install -y qemu \
			qemu binfmt-support qemu-user-static
		sudo docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
		DOCKER_IMAGE="${OS_ARCH}/${DOCKER_IMAGE}"
	fi
	echo 'DOCKER_OPTS="-H tcp://127.0.0.1:2375 -H unix:///var/run/docker.sock -s devicemapper"' | sudo tee /etc/default/docker > /dev/null
	sudo service docker restart
	sudo docker pull "$DOCKER_IMAGE"
}

__save_env_for_docker() {
	local env_file="$1/inside-travis-ci-docker-env"
	for env in $INSIDE_DOCKER_TRAVIS_CI_ENV ; do
		val="$(eval echo "\$${env}")"
		if [ -n "$val" ] ; then
			echo "export ${env}=\"${val}\"" >> "${env_file}"
		fi
	done
}

run_docker_script() {
	local DOCKER_SCRIPT="$(get_script_path "$1")"
	local MOUNTPOINT="${INSIDE_DOCKER_BUILD_DIR}"
	local DOCKER_IMAGE="${OS_TYPE}:${OS_VERSION}"

	if [ -n "$OS_ARCH" ] ; then
		DOCKER_IMAGE="${OS_ARCH}/${DOCKER_IMAGE}"
	fi

	__save_env_for_docker "$(pwd)"

	sudo docker run --rm=true \
		-v "$(pwd):/${MOUNTPOINT}:rw" \
		"$DOCKER_IMAGE" \
		/bin/bash -e "/${MOUNTPOINT}/${DOCKER_SCRIPT}" "${MOUNTPOINT}" "${OS_TYPE}"
}

command_exists() {
	local cmd=$1
	[ -n "$cmd" ] || return 1
	type "$cmd" >/dev/null 2>&1
}

ensure_command_exists() {
	local cmd="$1"
	local package="$2"
	local yes_confirm
	[ -n "$cmd" ] || return 1
	[ -n "$package" ] || package="$cmd"
	! command_exists "$cmd" || return 0
	# go through known package managers
	for pacman in apt-get brew yum ; do
		command_exists $pacman || continue
		if [ "$pacman" = "brew" ] ; then
			yes_confirm=
		else
			yes_confirm="-y"
		fi
		"$pacman" install $yes_confirm "$package" || {
			# Try an update if install doesn't work the first time
			"$pacman" $yes_confirm update && \
				"$pacman" install $yes_confirm "$package"
		}
		return $?
	done
	return 1
}

version_gt() { test "$(echo "$@" | tr " " "\n" | sort -V | head -n 1)" != "$1"; }
version_le() { test "$(echo "$@" | tr " " "\n" | sort -V | head -n 1)" = "$1"; }
version_lt() { test "$(echo "$@" | tr " " "\n" | sort -rV | head -n 1)" != "$1"; }
version_ge() { test "$(echo "$@" | tr " " "\n" | sort -rV | head -n 1)" = "$1"; }

get_codename() {
	local VERSION_CODENAME
	eval $(grep -w VERSION_CODENAME /etc/os-release)
	echo "$VERSION_CODENAME"
}

get_dist_id() {
	local ID
	eval $(grep -w ID /etc/os-release)
	echo "$ID"
}

get_version() {
	local VERSION_ID
	eval $(grep -w VERSION_ID /etc/os-release)
	echo "$VERSION_ID"
}

is_ubuntu_at_least_ver() {
	[ "$(get_dist_id)" = "ubuntu" ] || return 1
	version_ge "$(get_version)" "$1"
}

is_centos_at_least_ver() {
	[ "$(get_dist_id)" = "centos" ] || return 1
	version_ge "$(get_version)" "$1"
}

is_python_at_least_ver() {
	local out python_exec

	python_exec="$1"
	command_exists "$python_exec" || return 1
	out=$($python_exec --version)
	version_ge "${out#* }" "$2"
}

is_arm() {
	[ "$(dpkg --print-architecture)" = "armhf" ] || return 1
}

is_arm64() {
	[ "$(dpkg --print-architecture)" = "arm64" ] || return 1
}

print_github_api_rate_limits() {
	# See https://developer.github.com/v3/rate_limit/
	# Note: Accessing this endpoint does not count against your REST API rate limit.
	echo_green '-----------------------------------------'
	echo_green 'Github API Rate limits'
	echo_green '-----------------------------------------'
	wget -q -O- https://api.github.com/rate_limit
	echo_green '-----------------------------------------'
}

setup_build_type_env_vars() {
	OS_TYPE=${OS_TYPE:-default}

	# For a 'arm32_v7/debian_docker' string, OS TYPE becomes 'debian'
	# This also works for just 'debian_docker'
	# And we're extracting OS_ARCH if present
	if [ "${OS_TYPE#*_}" = "docker" ] ; then
		BUILD_TYPE=generic_docker
		OS_TYPE=${OS_TYPE%_*}
		OS_ARCH=${OS_TYPE%/*}
		OS_TYPE=${OS_TYPE#*/}
		if [ "$OS_ARCH" = "$OS_TYPE" ] ; then
			OS_ARCH=
		fi
	else
		BUILD_TYPE="$OS_TYPE"
	fi

	export OS_TYPE
	export OS_ARCH
	export BUILD_TYPE
}

ensure_command_exists wget
ensure_command_exists sudo

if [ -z "${LDIST}" -a -f "build/.LDIST" ] ; then
	export LDIST="-$(cat build/.LDIST)"
	echo $LDIST
fi
if [ -z "${LDIST}" ] || [ "$LDIST" = "DO_NOT_DEPLOY" ] ; then
	export LDIST="-$(get_ldist)"
        echo $LDIST
fi
