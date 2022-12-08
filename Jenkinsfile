#!groovy

pipeline {
  agent {
    label 'docker-runner'
  }
  stages{
    stage('Build docker image') {
      environment {
        AWS_ACCESS_KEY_ID = credentials('AWS_ACCESS_KEY_ID')
        AWS_SECRET_ACCESS_KEY = credentials('AWS_SECRET_ACCESS_KEY')
      }
      steps {
        echo "Building Docker image"
        sh """
        \$(aws ecr get-login --no-include-email --region us-west-2 &> /dev/null)
        docker build -t pkpdposterior .
        """
      }

    }
    stage('Test PKPDposterior') {
      steps {
        echo 'Installing and checking PKPDposterior'
        sh """
        docker run -d -t --name ${BUILD_TAG} pkpdposterior:latest
        docker exec -i ${BUILD_TAG} Rscript -e "devtools::check(vignettes = FALSE)"
        """
      }
    }
  }
  post {
    always {
      sh """
      docker rm -f ${BUILD_TAG} &>/dev/null && echo 'Removed container'
      """
    }
  }
}
