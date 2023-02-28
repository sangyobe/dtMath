
# GitLab CI/CD를 이용한 컴파일/빌드 자동화

CI/CD를 통해 개발 pc가 아닌 곳에서 컴파일과 빌드가 정상 동작하는지 확인하기 위해 아래의 문서를 작성합니다.



## Group Runner 이용 (등록과정 생략)
GitLab CI/CD를 이용하기 위해서는 Runner를 등록해야 하는데, Ubuntu의 경우 로보틱스랩에서 등록한 shared runner를 이용하는 것으로 등록 없이 바로 CI/CD 이용이 가능합니다. (개인적인 Runner 등록은 아래의 참조 사이트들을 보시면 됩니다.)
따라서, Ubuntu에서 개발자에게 필요한 최소의 컴파일/빌드 자동화 과정만 설명합니다.

아래와 같이, 자동화 하려는 GitLab 프로젝트에서 settings → CI/CD → Runners Expand → Group runners 에 등록된 ros-noetic 을 이용합니다.


![runner1](doc/image/cicd_runner_1.png)

![runner2](doc/image/cicd_runner_2.png)



## 프로젝트 폴더에서 '.gitlab-ci.yml' 파일 생성
GitLab은 runner 등록이 되어 있고, .gitlab-ci.yml 파일이 있으면  자동으로 CI/CD를 실행하기 때문에, .gitlab-ci.yml 을 작성하는 것으로 작업을 마무리할 수 있습니다. 여러 기능 없이 GitLab 프로젝트의 여러 commit들 중 tag가 업데이트 될때마다 compile/build를 수행하도록 작성하겠습니다. 예시 코드는 dtMath_cmake 프로젝트에서 확인 가능합니다.

a. 아래 파일 작성후 업로드
b. 원하는 commit에 버전 정보를 표시한 tag 생성 (ex: v1.0.0)

코드에서 ' - ./bin/<executable name>' 에 실행 파일 이름 작성후 주석 풀면 실행까지 테스트해볼 수 있습니다. (프로그램이 자동으로 종료가 되는 경우 해당. while 문 등으로 프로세스를 계속 소유하고 있는 프로그램의 경우 해당 X)

```yml
# in .gitlab-ci.yml
stages:
  - test
  - build


# Job name 'build'
build_linux:
  rules:
    - if: $CI_COMMIT_TAG =~ /.*v.*/
      when: always
    - when: never
    # - when: always # temparay set until windows build complete (230226)
# Runs at 'build' stage
  stage: build
  tags:
    # - dtMath_CMake_tag_windows
    - ros-noetic
  script:
    - apt update && DEBIAN_FRONTEND=noninteractive apt install -y git-core
    - rm -rf build  # Guarantee that tests and coverage are using latest src and tests
    - mkdir build
    - cd build
    - cmake .. 
    - cmake --build . --config Release
#    - ./bin/<executable name>


include:
  - template: Security/SAST.gitlab-ci.yml
  - template: Code-Quality.gitlab-ci.yml

code_quality:
  rules:
    - if: $CI_COMMIT_BRANCH == "main" 
    - if: $CI_PIPELINE_SOURCE == "merge_request_event"
  services:
  tags:
    - cq-sans-dind

sast:
  stage: test
```

## Tag 생성

Tag 생성은 업로드 된 commit에 tag를 지정합니다.

![tag1](doc/image/cicd_tag_1.png)

![tag2](doc/image/cicd_tag_2.png)

![tag3](doc/image/cicd_tag_3.png)







## CI/CD 확인

GitLab 프로젝트 → CI/CD → Pipelines → Status 에서 passed 를 통해 동작 완료를 확인할 수 있으며 오류가 있으면 failed가 됩니다. passed를 누르면 pipeline에서 진행된 내용들을 확인할 수 있으며, 세부 내용 선택 시, gitlab-runner가 수행한 상세 작업과정들을 확인할 수 있고, 오류 발생시 에러메시지 또한 확인 가능합니다.


![pipeline1](doc/image/cicd_pipeline_1.png)

![pipeline2](doc/image/cicd_pipeline_2.png)

![pipeline3](doc/image/cicd_pipeline_3.png)











